import logging

import libsbml

from mod_sbml.sbml.ubiquitous_manager import UBIQUITOUS_THRESHOLD, select_metabolite_ids_by_term_ids, \
    get_ubiquitous_chebi_ids
from mod_sbml.onto import filter_ontology
from mod_sbml.sbml.compartment.compartment_manager import separate_boundary_metabolites
from mod_sbml.sbml.submodel_manager import get_biomass_r_ids
from sbml_generalization.sbml.sbml_helper import save_as_comp_generalized_sbml, remove_is_a_reactions, \
    remove_unused_elements
from sbml_generalization.generalization.model_generalizer import generalize_species, generalize_reactions
from mod_sbml.annotation.chebi.chebi_annotator import add_equivalent_chebi_ids, \
    EQUIVALENT_RELATIONSHIPS, annotate_metabolites, get_species_id2chebi_id
from mod_sbml.utils.misc import invert_map

__author__ = 'anna'


def get_ub_elements(model, chebi, s_id2chebi_id, ub_chebi_ids, ub_s_ids):
    """
    Infers ubiquitous species in the model.
    :param s_id2chebi_id: dict {species_id: ChEBI_term_id}
    :param model: libsbml.Model, input model
    :param chebi: mod_sbml.onto.obo_ontology.Ontology ChEBI ontology
    :param ub_s_ids: optional, ids of ubiquitous species (will be inferred if set to None)
    :param ub_chebi_ids: optional, ids of ubiquitous ChEBI terms (will be inferred if set to None)
    :return: tuple (ub_chebi_ids, ub_s_ids): set of ubiquitous ChEBI term ids, set of ubiquitous species ids.
    """
    if ub_s_ids:
        if not ub_chebi_ids:
            ub_chebi_ids = set()
        ub_chebi_ids |= {s_id2chebi_id[s_id] for s_id in ub_s_ids if s_id in s_id2chebi_id}
    else:
        if not ub_chebi_ids:
            ub_chebi_ids = get_ubiquitous_chebi_ids(add_common=True, add_cofactors=True, add_frequent=False,
                                                    chebi=chebi)
        else:
            ub_chebi_ids = add_equivalent_chebi_ids(chebi, ub_chebi_ids)
        ub_s_ids = select_metabolite_ids_by_term_ids(model, ub_chebi_ids)
    return ub_chebi_ids, ub_s_ids


def generalize_model(in_sbml, chebi, groups_sbml, out_sbml, ub_s_ids=None, ub_chebi_ids=None, ignore_biomass=True):
    """
    Generalizes a model.
    :param in_sbml: str, path to the input SBML file
    :param chebi: mod_sbml.onto.obo_ontology.Ontology ChEBI ontology
    :param groups_sbml: str, path to the output SBML file (with groups extension)
    :param out_sbml: str, path to the output SBML file (generalized)
    :param ub_s_ids: optional, ids of ubiquitous species (will be inferred if set to None)
    :param ub_chebi_ids: optional, ids of ubiquitous ChEBI terms (will be inferred if set to None)
    :param ignore_biomass: boolean, whether to ignore the biomass reaction (and its stoichiometry preserving constraint)
    :return: tuple (r_id2g_eq, s_id2gr_id, s_id2chebi_id, ub_s_ids):
    dict {reaction_id: reaction_group_id}, dict {species_id: species_group_id}, dict {species_id: ChEBI_term_id},
    collection of ubiquitous species_ids.
    """
    # input_model
    input_doc = libsbml.SBMLReader().readSBML(in_sbml)
    input_model = input_doc.getModel()
    r_ids_to_ignore = get_biomass_r_ids(input_model) if ignore_biomass else None

    remove_is_a_reactions(input_model)
    annotate_metabolites(input_model, chebi)
    # TODO: fix comp separation
    # separate_boundary_metabolites(input_model)
    remove_unused_elements(input_model)

    logging.info("mapping species to ChEBI")
    s_id2chebi_id = get_species_id2chebi_id(input_model)
    ub_chebi_ids, ub_s_ids = get_ub_elements(input_model, chebi, s_id2chebi_id, ub_chebi_ids, ub_s_ids)

    terms = (t for t in (chebi.get_term(t_id) for t_id in s_id2chebi_id.values()) if t)
    old_onto_len = len(chebi)
    filter_ontology(chebi, terms, relationships=EQUIVALENT_RELATIONSHIPS, min_deepness=3)
    logging.info('Filtered the ontology from %d terms to %d' % (old_onto_len, len(chebi)))

    threshold = min(max(3, int(0.1 * input_model.getNumReactions())), UBIQUITOUS_THRESHOLD)
    s_id2clu, ub_s_ids = generalize_species(input_model, s_id2chebi_id, ub_s_ids, chebi, ub_chebi_ids, threshold,
                                            r_ids_to_ignore=r_ids_to_ignore)
    logging.info("generalized species")
    r_id2clu = generalize_reactions(input_model, s_id2clu, s_id2chebi_id, ub_chebi_ids,
                                    r_ids_to_ignore=r_ids_to_ignore)
    logging.info("generalized reactions")

    clu2s_ids = {(c_id, term): s_ids for ((c_id, (term, )), s_ids) in invert_map(s_id2clu).items()}
    r_id2g_eq, s_id2gr_id = save_as_comp_generalized_sbml(input_model, out_sbml, groups_sbml, r_id2clu, clu2s_ids,
                                                          ub_s_ids, chebi)
    return r_id2g_eq, s_id2gr_id, s_id2chebi_id, ub_s_ids


def ubiquitize_model(in_sbml, chebi, groups_sbml, ub_s_ids=None, ub_chebi_ids=None):
    """
    Infers and marks ubiquitous species in the model.
    :param in_sbml: str, path to the input SBML file
    :param chebi: mod_sbml.onto.obo_ontology.Ontology ChEBI ontology
    :param groups_sbml: str, path to the output SBML file (with groups extension)
    :param ub_s_ids: optional, ids of ubiquitous species (will be inferred if set to None)
    :param ub_chebi_ids: optional, ids of ubiquitous ChEBI terms (will be inferred if set to None)
    :return: tuple (s_id2chebi_id, ub_s_ids): dict {species_id: ChEBI_term_id},  collection of ubiquitous species_ids.
    """
    input_doc = libsbml.SBMLReader().readSBML(in_sbml)
    input_model = input_doc.getModel()
    annotate_metabolites(input_model, chebi)

    logging.info("mapping species to ChEBI")
    s_id2chebi_id = get_species_id2chebi_id(input_model)
    _, ub_s_ids = get_ub_elements(input_model, chebi, s_id2chebi_id, ub_chebi_ids, ub_s_ids)

    save_as_comp_generalized_sbml(input_model, None, groups_sbml, {}, {}, ub_s_ids, chebi)
    return s_id2chebi_id, ub_s_ids
