import logging

import libsbml

from mod_sbml.annotation.chebi.chebi_annotator import get_chebi_term_by_annotation, CHEBI_PREFIX
from mod_sbml.utils.misc import invert_map
from mod_sbml.annotation.rdf_annotation_helper import add_annotation, get_annotations
from mod_sbml.sbml.sbml_manager import get_products, get_reactants, get_metabolites, generate_unique_id, create_reaction, \
    create_species

GROUP_TYPE_EQUIV = "equivalent"

GROUP_TYPE_UBIQUITOUS = "ubiquitous"

SBO_CHEMICAL_MACROMOLECULE = "SBO:0000248"

SBO_BIOCHEMICAL_REACTION = "SBO:0000176"


__author__ = 'anna'


def check(value, message):
    """
    If 'value' is None, logs an error message constructed using 'message'.
    If 'value' is an integer, it assumes it is a libSBML return status code.
    If the code value is LIBSBML_OPERATION_SUCCESS, returns without further action;
    if it is not, logs an error message constructed using 'message' along with text from
    libSBML explaining the meaning of the code.
    """
    if value is None:
        logging.error('LibSBML returned a null value trying to %s.' % message)
    elif type(value) is int and value != libsbml.LIBSBML_OPERATION_SUCCESS:
        logging.error('Error encountered trying to %s. LibSBML returned error code %s: "%s".'
                      % (message, value, libsbml.OperationReturnValue_toString(value).strip()))


def create_model():
    """Returns a document with an empty SBML Level 3 model."""

    try:
        document = libsbml.SBMLDocument(3, 1)
    except ValueError:
        raise SystemExit('Could not create SBMLDocument object')

    # Create the basic Model object inside the SBMLDocument object.  To
    # produce a model with complete units for the reaction rates, we need
    # to set the 'timeUnits' and 'extentUnits' attributes on Model.  We
    # set 'substanceUnits' too, for good measure, though it's not strictly
    # necessary here because we also set the units for invididual species
    # in their definitions.
    model = document.createModel()
    check(model, 'create model')
    check(model.setTimeUnits("second"), 'set model-wide time units')
    check(model.setExtentUnits("mole"), 'set model units of extent')
    check(model.setSubstanceUnits('mole'), 'set model substance units')

    # Create a unit definition we will need later.  Note that SBML Unit
    # objects must have all four attributes 'kind', 'exponent', 'scale'
    # and 'multiplier' defined.
    per_second = model.createUnitDefinition()
    check(per_second, 'create unit definition')
    check(per_second.setId('per_second'), 'set unit definition id')
    unit = per_second.createUnit()
    check(unit, 'create unit on per_second')
    check(unit.setKind(libsbml.UNIT_KIND_SECOND), 'set unit kind')
    check(unit.setExponent(-1), 'set unit exponent')
    check(unit.setScale(0), 'set unit scale')
    check(unit.setMultiplier(1), 'set unit multiplier')

    return document


def remove_unused_elements(model, include_modifiers=True):
    """
    Removes unused element from the model, i.e. compartments that contain no metabolites,
    and metabolites that do not participate in any reaction.
    :param include_modifiers: whether to keep reaction modifiers.
    :param model: libsbml.Model model of interest
    :return: void, input model is modified inplace
    """
    species_to_keep = get_used_species(model, include_modifiers)
    sp_list = list(model.getListOfSpecies())
    for species in sp_list:
        species_id = species.getId()
        if not (species_id in species_to_keep):
            model.removeSpecies(species_id)
    compartments_to_keep = set()
    for species in model.getListOfSpecies():
        compartment_id = species.getCompartment()
        compartments_to_keep.add(compartment_id)
        comp = model.getCompartment(compartment_id)
        if comp:
            outer_compartment = comp.getOutside()
            while outer_compartment:
                compartments_to_keep.add(outer_compartment)
                outer_compartment = model.getCompartment(outer_compartment).getOutside()
    c_list = list(model.getListOfCompartments())
    for compartment in c_list:
        compartment_id = compartment.getId()
        if not (compartment_id in compartments_to_keep):
            model.removeCompartment(compartment_id)


def get_used_species(model, include_modifiers=True):
    """
    Returns a list of species participating in any of the model reactions.
    :param include_modifiers: whether to include reaction modifiers
    :param model: libsbml.Model model of interest
    :return: set of used species ids
    """
    species_to_keep = set()
    for reaction in model.getListOfReactions():
        species_to_keep |= get_metabolites(reaction, include_modifiers=include_modifiers)
    return species_to_keep


def normalize(t):
    if isinstance(t, tuple):
        t = ''.join(flatten(t))
    return ''.join(e for e in t if e.isalnum()).lower()


def flatten(t):
    if not isinstance(t, tuple):
        return t,
    elif not len(t):
        return ()
    else:
        return flatten(t[0]) + flatten(t[1:])


def remove_is_a_reactions(model):
    """
    Removes 'isa' reactions from the model. 'isa' reactions can be found in yeast.net models and contain
    purely hierarchical information, e.g. 'lauroyl-CoA isa fatty-acyl-CoA'.
    :param model: libsbml.Model model of interest
    :return: void (input model is modified inplace)
    """
    to_remove = []
    for reaction in model.getListOfReactions():
        if 1 == reaction.getNumReactants() == reaction.getNumProducts() \
                and reaction.getName().find("isa ") != -1 \
                and model.getCompartment(reaction.getListOfReactants().get(0).getSpecies()) == \
                        model.getCompartment(reaction.getListOfReactants().get(0).getSpecies()):
            to_remove.append(reaction.getId())
    for r_id in to_remove:
        model.removeReaction(r_id)


def set_consistency_level(doc):
    for consistency in [libsbml.LIBSBML_CAT_GENERAL_CONSISTENCY,
                        libsbml.LIBSBML_CAT_IDENTIFIER_CONSISTENCY,
                        libsbml.LIBSBML_CAT_UNITS_CONSISTENCY,
                        libsbml.LIBSBML_CAT_MATHML_CONSISTENCY,
                        libsbml.LIBSBML_CAT_SBO_CONSISTENCY,
                        libsbml.LIBSBML_CAT_OVERDETERMINED_MODEL,
                        libsbml.LIBSBML_CAT_MODELING_PRACTICE]:
        doc.setConsistencyChecksForConversion(consistency, False)


def convert_to_lev3_v1(model):
    doc = libsbml.SBMLDocument(model.getSBMLNamespaces())
    doc.setModel(model)
    set_consistency_level(doc)
    doc.checkL3v1Compatibility()
    doc.setLevelAndVersion(3, 1, False)
    doc.enablePackage(libsbml.GroupsExtension.getXmlnsL3V1V1(), "groups", True)
    doc.enablePackage(libsbml.LayoutExtension.getXmlnsL3V1V1(), "layout", True)
    # doc.getSBMLNamespaces().addPackageNamespace("groups", 1)
    return doc


def save_as_comp_generalized_sbml(input_model, out_sbml, groups_sbml, r_id2clu, clu2s_ids, ub_sps, onto):
    logging.info("serializing generalization")
    s_id_increment, r_id_increment = 0, 0

    if groups_sbml:
        doc = convert_to_lev3_v1(input_model)
        groups_model = doc.getModel()
        groups_plugin = groups_model.getPlugin("groups")
        if groups_plugin:
            logging.info("  saving ubiquitous species annotations")
            s_group = groups_plugin.createGroup()
            s_group.setId("g_ubiquitous_sps")
            s_group.setKind(libsbml.GROUP_KIND_COLLECTION)
            s_group.setSBOTerm(SBO_CHEMICAL_MACROMOLECULE)
            s_group.setName("ubiquitous species")
            for s_id in ub_sps:
                member = s_group.createMember()
                member.setIdRef(s_id)
            add_annotation(s_group, libsbml.BQB_IS_DESCRIBED_BY, GROUP_TYPE_UBIQUITOUS)
    if out_sbml:
        # generalized model
        generalized_doc = convert_to_lev3_v1(input_model)
        generalized_model = generalized_doc.getModel()
        for _ in range(0, generalized_model.getNumReactions()):
            generalized_model.removeReaction(0)

    r_id2g_eq, s_id2gr_id = {}, {}
    if not clu2s_ids:
        logging.info("  nothing to serialize")
    else:
        clu2r_ids = invert_map(r_id2clu)
        logging.info("  creating species groups")
        for ((c_id, t), s_ids) in clu2s_ids.items():
            comp = input_model.getCompartment(c_id)
            if len(s_ids) > 1:
                t = onto.get_term(t)
                t_name, t_id = (t.get_name(), t.get_id()) if t \
                    else (' or '.join(input_model.getSpecies(s_id).getName() for s_id in s_ids), None)
                if not t_id:
                    t = t_name

                if out_sbml:
                    new_species = create_species(model=generalized_model, compartment_id=comp.getId(), type_id=None,
                                                 name="{0} ({1}) [{2}]".format(t_name, len(s_ids), comp.getName()))
                    add_annotation(new_species, libsbml.BQB_IS, t_id, CHEBI_PREFIX)
                    new_s_id = new_species.getId()
                else:
                    s_id_increment += 1
                    new_s_id = generate_unique_id(input_model, "s_g_", s_id_increment)
                for s_id in s_ids:
                    s_id2gr_id[s_id] = new_s_id, t

                if groups_sbml and groups_plugin:
                    # save as a group
                    s_group = groups_plugin.createGroup()
                    s_group.setId(new_s_id)
                    s_group.setKind(libsbml.GROUP_KIND_CLASSIFICATION)
                    s_group.setSBOTerm(SBO_CHEMICAL_MACROMOLECULE)
                    g_name = "{0} [{1}]".format(t_name, comp.getName())
                    s_group.setName(g_name)
                    # logging.info("%s: %d" % (g_name, len(s_ids)))
                    if t_id:
                        add_annotation(s_group, libsbml.BQB_IS, t_id, CHEBI_PREFIX)
                    for s_id in s_ids:
                        member = s_group.createMember()
                        member.setIdRef(s_id)
                    add_annotation(s_group, libsbml.BQB_IS_DESCRIBED_BY, GROUP_TYPE_EQUIV)

        generalize_species = lambda species_id: s_id2gr_id[species_id][0] if (species_id in s_id2gr_id) else species_id
        s_id_to_generalize = set(s_id2gr_id.keys())
        logging.info("  creating reaction groups")
        for clu, r_ids in clu2r_ids.items():
            representative = input_model.getReaction(list(r_ids)[0])
            r_name = "generalized %s" % representative.getName()
            if out_sbml:
                reactants = dict(get_reactants(representative, stoichiometry=True))
                products = dict(get_products(representative, stoichiometry=True))
                if (len(r_ids) == 1) and \
                        not ((set(reactants.keys()) | set(products.keys())) & s_id_to_generalize):
                    create_reaction(generalized_model, reactants, products, name=representative.getName(),
                                    reversible=representative.getReversible(), id_=representative.getId())
                    continue
                r_id2st = {generalize_species(it): st for (it, st) in reactants.items()}
                p_id2st = {generalize_species(it): st for (it, st) in products.items()}
                reversible = next((False for r_id in r_ids if not input_model.getReaction(r_id).getReversible()), True)
                new_r_id = create_reaction(generalized_model, r_id2st, p_id2st, name=r_name, reversible=reversible,
                                           id_=representative.getId() if len(r_ids) == 1 else None).getId()
            elif len(r_ids) > 1:
                r_id_increment += 1
                new_r_id = generate_unique_id(input_model, "r_g_", r_id_increment)
            if len(r_ids) > 1:
                for r_id in r_ids:
                    r_id2g_eq[r_id] = new_r_id, r_name
                if groups_sbml and groups_plugin:
                    # save as a group
                    r_group = groups_plugin.createGroup()
                    r_group.setId(new_r_id)
                    r_group.setKind(libsbml.GROUP_KIND_COLLECTION)
                    r_group.setSBOTerm(SBO_BIOCHEMICAL_REACTION)
                    r_group.setName(r_name)
                    for r_id in r_ids:
                        member = r_group.createMember()
                        member.setIdRef(r_id)
                    add_annotation(r_group, libsbml.BQB_IS_DESCRIBED_BY, GROUP_TYPE_EQUIV)
    if out_sbml:
        remove_unused_elements(generalized_model)
        save_as_sbml(generalized_model, out_sbml)
    if groups_sbml and groups_model:
        save_as_sbml(groups_model, groups_sbml)

    logging.info("serialized to " + groups_sbml)
    return r_id2g_eq, s_id2gr_id


def save_as_sbml(input_model, out_sbml):
    logging.info("saving to {0}".format(out_sbml))
    out_doc = libsbml.SBMLDocument(input_model.getSBMLNamespaces())
    out_doc.setModel(input_model)
    libsbml.writeSBMLToFile(out_doc, out_sbml)


def parse_group_sbml(groups_sbml, chebi):
    doc = libsbml.SBMLReader().readSBMLFromFile(groups_sbml)
    groups_model = doc.getModel()
    groups_plugin = groups_model.getPlugin("groups")
    r_id2g_id, s_id2gr_id, ub_sps = {}, {}, set()
    if groups_plugin:
        for group in groups_plugin.getListOfGroups():
            gr_members = [it.getIdRef() for it in group.getListOfMembers()]
            gr_id, gr_name = group.getId(), group.getName()
            gr_sbo = group.getSBOTermID()
            gr_type = next(get_annotations(group, libsbml.BQB_IS_DESCRIBED_BY), None)
            if not gr_type:
                continue
            if SBO_BIOCHEMICAL_REACTION == gr_sbo:
                if GROUP_TYPE_EQUIV == gr_type:
                    for r_id in gr_members:
                        r_id2g_id[r_id] = gr_id, gr_name, len(gr_members)
            elif SBO_CHEMICAL_MACROMOLECULE == gr_sbo:
                if GROUP_TYPE_UBIQUITOUS == gr_type:
                    ub_sps = set(gr_members)
                elif GROUP_TYPE_EQUIV == gr_type:
                    for s_id in gr_members:
                        term = get_chebi_term_by_annotation(group, chebi)
                        s_id2gr_id[s_id] = gr_id, term if term else gr_name, len(gr_members)
    else:
        raise GrPlError()
    return r_id2g_id, s_id2gr_id, ub_sps


def check_for_groups(groups_sbml, sbo_term, group_type):
    doc = libsbml.SBMLReader().readSBMLFromFile(groups_sbml)
    groups_plugin = doc.getModel().getPlugin("groups")
    if groups_plugin:
        for group in groups_plugin.getListOfGroups():
            gr_sbo = group.getSBOTermID()
            gr_type = next(get_annotations(group, libsbml.BQB_IS_DESCRIBED_BY), None)
            if sbo_term == gr_sbo and group_type == gr_type:
                return True
    return False


class GrPlError(Exception):
    def __init__(self):
        self.msg = "groups plugin not installed"

