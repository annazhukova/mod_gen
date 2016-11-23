from collections import defaultdict
from itertools import chain
from mod_sbml.sbml.sbml_manager import get_products, get_reactants
from mod_sbml.sbml.ubiquitous_manager import get_proton_ch_ids

__author__ = 'anna'


def get_vertical_key(model, r, s_id2clu, s_id2term_id, ubiquitous_chebi_ids):
    """
    Gets a reaction key: ubiquitous_reactants, ubiquitous_products,
    specific_reactant_classes, specific_product_classes
    (for reversible reactions reactants and products can be swapped according to sorting)
    :param model: libsbml.Model
    :param r: libsbml.Reaction reaction of interest
    :param s_id2clu: dict {metabolite_id: (compartment_id, cluster)}
    :param s_id2term_id: dict {metabolite_id: ChEBI_term_id}
    :param ubiquitous_chebi_ids: set of ubiquitous ChEBI_ids
    :return: tuple (ubiquitous_reactants, ubiquitous_products,
    specific_reactant_classes, specific_product_classes)
    """
    ubiquitous_reactants, ubiquitous_products, specific_reactant_classes, specific_product_classes = \
        get_key_elements(model, r, s_id2clu, s_id2term_id, ubiquitous_chebi_ids)
    if r.getReversible() and need_to_reverse(
            ubiquitous_reactants, ubiquitous_products, specific_reactant_classes, specific_product_classes):
        return ubiquitous_products, ubiquitous_reactants, specific_product_classes, specific_reactant_classes
    return ubiquitous_reactants, ubiquitous_products, specific_reactant_classes, specific_product_classes


def need_to_reverse(ubiquitous_reactants, ubiquitous_products, specific_reactant_classes, specific_product_classes):
    return (ubiquitous_reactants > ubiquitous_products) or (
        not ubiquitous_reactants and not ubiquitous_products and (
            len(specific_reactant_classes) > len(specific_product_classes) or (
                len(specific_reactant_classes) == len(
                    specific_product_classes) and specific_reactant_classes > specific_product_classes)))


def get_r_compartments(model, r):
    return tuple({model.getSpecies(s_id).getCompartment()
                  for s_id in chain((species_ref.getSpecies() for species_ref in r.getListOfReactants()),
                                    (species_ref.getSpecies() for species_ref in r.getListOfProducts()))})


def vertical_key2simplified_vertical_key(vertical_key):
    """
    Simplifies a given reaction key by replacing specific reaction and product clusters by their numbers.
    :param vertical_key: tuple (ubiquitous_reactants, ubiquitous_products,
    specific_reactant_classes, specific_product_classes) -- reaction key
    :return: simplified vertical key (ubiquitous_reactants, ubiquitous_products,
    number_of_specific_reactant_classes, number_of_specific_product_classes)
    """
    u_rs, u_ps, rs, ps = vertical_key
    s_rs = u_rs, len(rs)
    s_ps = u_ps, len(ps)
    if s_rs > s_ps:
        s_rs, s_ps = s_ps, s_rs
    return s_rs, s_ps


def get_vk2r_ids(model, s_id2clu, s_id2term_id, ubiquitous_chebi_ids, r_ids_to_ignore=None):
    """
    Calculates key to reaction ids mapping based on the metabolite clustering.
    :param model: libsbml.Model model of interest
    :param s_id2clu: dict {metabolite_id: (compartment_id, cluster)}
    :param s_id2term_id: dict {metabolite_id: ChEBI_term_id}
    :param ubiquitous_chebi_ids: set of ubiquitous ChEBI_ids
    :param r_ids_to_ignore: (optional) ids of reactions whose stoichiometry preserving constraint can be ignores
    :return: dict {key: reaction_id_set}
    """
    vk2r = defaultdict(set)
    for r in model.getListOfReactions():
        if r_ids_to_ignore and r.getId() in r_ids_to_ignore:
            continue
        vk2r[get_vertical_key(model, r, s_id2clu, s_id2term_id, ubiquitous_chebi_ids)].add(r.getId())
    return vk2r


def is_reactant(model, t_id, r, s_id2clu, s_id2term_id, ubiquitous_chebi_ids):
    ubiquitous_reactants, ubiquitous_products, specific_reactant_classes, specific_product_classes = \
        get_key_elements(model, r, s_id2clu, s_id2term_id, ubiquitous_chebi_ids)
    if r.getReversible() and need_to_reverse(
            ubiquitous_reactants, ubiquitous_products, specific_reactant_classes, specific_product_classes):
        return t_id in {s_id2term_id[s_id] if s_id in s_id2term_id else s_id for s_id in get_products(r)}
    else:
        return t_id in {s_id2term_id[s_id] if s_id in s_id2term_id else s_id for s_id in get_reactants(r)}


def get_key_elements(model, r, s_id2clu, s_id2term_id, ubiquitous_chebi_ids, ignore_ch_ids=get_proton_ch_ids()):
    """
    Gets elements that compose a reaction key: ubiquitous_reactants, ubiquitous_products,
    specific_reactant_classes, specific_product_classes
    :param model: libsbml.Model
    :param r: libsbml.Reaction reaction of interest
    :param s_id2clu: dict {metabolite_id: (compartment_id, cluster)}
    :param s_id2term_id: dict {metabolite_id: ChEBI_term_id}
    :param ubiquitous_chebi_ids: set of ubiquitous ChEBI_ids
    :param ignore_ch_ids: set of ChEBI_ids to be excluded from the result (by default protons)
    if there is anything else in the result
    :return: tuple (ubiquitous_reactants, ubiquitous_products,
    specific_reactant_classes, specific_product_classes)
    """

    def classify(s_ids):
        specific, ubiquitous, ignored_ubs = [], [], []
        for s_id in s_ids:
            c_id = model.getSpecies(s_id).getCompartment()
            if ubiquitous_chebi_ids and s_id in s_id2term_id and s_id2term_id[s_id] in ubiquitous_chebi_ids:
                if s_id2term_id[s_id] in ignore_ch_ids:
                    ignored_ubs.append((s_id2term_id[s_id], c_id))
                else:
                    ubiquitous.append((s_id2term_id[s_id], c_id))
            else:
                specific.append((s_id2clu[s_id][1] if s_id in s_id2clu
                                else ((s_id2term_id[s_id] if s_id in s_id2term_id else s_id), ), c_id))
        transform = lambda collection: tuple(sorted(collection))
        return transform(specific), transform(ubiquitous), transform(ignored_ubs)

    specific_reactant_classes, ubiquitous_reactants, ignored_reactants = classify(get_reactants(r))
    specific_product_classes, ubiquitous_products, ignored_products = classify(get_products(r))
    if not ubiquitous_reactants and not ubiquitous_products \
            and not specific_reactant_classes and not specific_product_classes:
        ubiquitous_reactants, ubiquitous_products = ignored_reactants, ignored_products
    return ubiquitous_reactants, ubiquitous_products, specific_reactant_classes, specific_product_classes

