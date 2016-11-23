from collections import Counter
from itertools import chain
import logging
from mod_sbml.annotation.chebi.chebi_annotator import EQUIVALENT_RELATIONSHIPS

from sbml_generalization.generalization.MaximizingThread import MaximizingThread
from sbml_generalization.generalization.StoichiometryFixingThread import StoichiometryFixingThread, compute_s_id2clu, \
    infer_clusters, suggest_clusters
from sbml_generalization.generalization.vertical_key import get_vk2r_ids
from mod_sbml.utils.misc import invert_map
from mod_sbml.onto.term import Term
from mod_sbml.sbml.ubiquitous_manager import UBIQUITOUS_THRESHOLD, get_frequent_term_ids, \
    select_metabolite_ids_by_term_ids

__author__ = 'anna'


def generalize_reactions(model, s_id2clu, s_id2term_id, ubiquitous_chebi_ids, r_ids_to_ignore=None):
    vk2r = get_vk2r_ids(model, s_id2clu, s_id2term_id, ubiquitous_chebi_ids, r_ids_to_ignore=r_ids_to_ignore)
    r_id2clu, i = {}, 0
    for r_ids in vk2r.values():
        for r_id in r_ids:
            r_id2clu[r_id] = i
        i += 1
    return r_id2clu


def maximize(unmapped_s_ids, model, term_id2clu, species_id2term_id, ub_chebi_ids, r_ids_to_ignore=None):
    clu2term_ids = invert_map(term_id2clu)
    s_id2clu = compute_s_id2clu(unmapped_s_ids, model, species_id2term_id, term_id2clu)

    r_id2clu = generalize_reactions(model, s_id2clu, species_id2term_id, ub_chebi_ids, r_ids_to_ignore=r_ids_to_ignore)

    thrds = []
    for (clu, term_ids) in clu2term_ids.items():
        if len(term_ids) <= 1:
            continue

        thread = MaximizingThread(model, term_ids, species_id2term_id, clu, term_id2clu,
                                  s_id2clu, ub_chebi_ids, r_id2clu, r_ids_to_ignore=r_ids_to_ignore)
        thrds.append(thread)
        thread.start()  # This actually causes the thread to run
    for th in thrds:
        th.join()  # This waits until the thread has completed
    return term_id2clu


def cover_t_ids(model, species_id2term_id, ubiquitous_t_ids, t_ids, onto, clu=None, r_ids_to_ignore=None):
    """
    Find ancestor terms that cover (generalize) given terms.
    :param model: libsbml.Model model of interest
    :param species_id2term_id: dict {species_id: term_id}
    :param ubiquitous_t_ids: collection of ubiquitous term ids
    :param t_ids: collection of term ids to be covered
    :param onto: mod_sbml.onto.obo_ontology.Ontology ChEBI ontology
    :param clu: current cluster to which the terms belong
    :param r_ids_to_ignore: collection of reaction ids to ignore (don't fix their Stoichiometry preserving constraints)
    :return: dictionary {term_id: cluster}
    """
    term_id2clu = {}
    real_terms = {onto.get_term(t_id) for t_id in t_ids if onto.get_term(t_id)}

    # If there is no term for t_id in the ontology, we assume it is a metabolite id instead
    unmapped_s_ids = {t_id for t_id in t_ids if not onto.get_term(t_id)}

    roots = onto.common_points(real_terms, relationships=EQUIVALENT_RELATIONSHIPS)
    if roots:
        root_id = roots[0].get_id()
        new_clu = clu + (root_id, ) if clu else (root_id, )
        return {t_id: new_clu for t_id in t_ids}

    roots = set()
    for term in real_terms:
        roots |= onto.get_generalized_ancestors_of_level(term, set(), None, 4)
    terms2root = {tuple(sorted(t.get_id() for t in onto.get_sub_tree(root))): root.get_id() for root in roots}
    for t_set, root_id in greedy({t.get_id() for t in real_terms}, terms2root, {it: 1 for it in terms2root}):
        new_clu = clu + (root_id, ) if clu else (root_id, )
        term_id2clu.update({t_id: new_clu for t_id in t_set})

    s_id2clu = compute_s_id2clu(set(), model, species_id2term_id, term_id2clu)
    infer_clusters(model, unmapped_s_ids, s_id2clu, species_id2term_id, ubiquitous_t_ids,
                   r_ids_to_ignore=r_ids_to_ignore)
    # for s_id in unmapped_s_ids:
    #     if s_id in s_id2clu:
    #         term_id2clu[s_id] = (s_id2clu[s_id][1], )
    return term_id2clu


def update_onto(onto, term_id2clu):
    ancestors = []
    clu2t_ids = invert_map(term_id2clu)
    for clu, t_ids in clu2t_ids.items():
        if len(t_ids) <= 1:
            continue
        terms = {onto.get_term(t_id) for t_id in t_ids if onto.get_term(t_id)}
        if terms:
            ancestors.extend(set(onto.common_points(terms, relationships=EQUIVALENT_RELATIONSHIPS)))
    removed_something = False
    count = Counter(ancestors)
    for t in (t for t in count.keys() if count[t] > 1):
        # if this term has been already removed as an ancestor/equivalent of another term
        if not onto.get_term(t.get_id()):
            continue
        for it in onto.get_generalized_ancestors(t, relationships=EQUIVALENT_RELATIONSHIPS):
            onto.remove_term(it, True)
        for it in onto.get_equivalents(t, relationships=EQUIVALENT_RELATIONSHIPS):
            onto.remove_term(it, True)
        onto.remove_term(t, True)
        removed_something = True
    return removed_something


def fix_stoichiometry(model, term_id2clu, species_id2term_id, ub_chebi_ids, onto, r_ids_to_ignore=None):
    clu2term_ids = invert_map(term_id2clu)
    thrds = []
    conflicts = []
    for r in model.getListOfReactions():
        if r_ids_to_ignore and r.getId() in r_ids_to_ignore:
            continue
        t_ids = {species_id2term_id[s_id] if s_id in species_id2term_id else s_id
                 for s_id in chain((species_ref.getSpecies() for species_ref in r.getListOfReactants()),
                                   (species_ref.getSpecies() for species_ref in r.getListOfProducts()))}
        if len(t_ids) > 1:
            conflicts.append(t_ids)
    for clu, term_ids in clu2term_ids.items():
        if len(term_ids) <= 1:
            continue
        clu_conflicts = [set(it) for it in {tuple(t_ids & term_ids) for t_ids in conflicts} if len(it) > 1]
        real_term_ids = {t_id for t_id in term_ids if onto.get_term(t_id)}
        unmapped_s_ids = {s_id for s_id in term_ids if not onto.get_term(s_id)}
        if clu_conflicts:
            thread = StoichiometryFixingThread(model, species_id2term_id, ub_chebi_ids, unmapped_s_ids, real_term_ids,
                                               clu_conflicts, onto, clu, term_id2clu, r_ids_to_ignore=r_ids_to_ignore)
            thrds.append(thread)
            thread.start()  # This actually causes the thread to run
    for th in thrds:
        th.join()  # This waits until the thread has completed


def greedy(yet_to_be_covered, set2label, set2score):
    """
    Greedy set coverage.
    :param yet_to_be_covered: set of interest
    :param set2label: dict {set: label} available sets and their labels
    :param set2score: dict {set: score}
    :return: iterator of tuples (set, label) corresponding to the selected sets for the coverage
    """
    yet_to_be_covered = set(yet_to_be_covered)
    while yet_to_be_covered and set2label:
        s = max(set2label.keys(),
                key=lambda candidate_terms: (len(set(candidate_terms) & yet_to_be_covered), set2score[candidate_terms]))
        result = set(s)
        # yield result
        yield result & yet_to_be_covered, set2label[s]
        yet_to_be_covered -= result
        del set2label[s]


def select_representative_terms(term_id2clu, onto):
    """
    Replaces each cluster with a tuple containing one ChEBI term id: (term_id, )
    :param term_id2clu: dict {term_id: cluster}
    :param onto: mod_sbml.onto.obo_ontology.Ontology ChEBI ontology
    :return: updated (inplace) dict term_id2clu {term_id: (cluster_representative_term_id, )}
    """
    clu2term_ids = invert_map(term_id2clu)
    used = set()
    i = 0
    for clu, term_ids in clu2term_ids.items():
        terms = {onto.get_term(t) for t in term_ids if onto.get_term(t)}
        common_ancestors = \
            {t for t in onto.common_points(terms, relationships=EQUIVALENT_RELATIONSHIPS)} if terms else set()
        options = common_ancestors - used
        if options:
            common_ancestor_term = options.pop()
        else:
            name = common_ancestors.pop().get_name() + " (another)" if common_ancestors else 'fake term'
            common_ancestor_term = Term(onto=onto, t_id="chebi:unknown_{0}".format(i), name=name)
            onto.add_term(common_ancestor_term)
            i += 1
        used.add(common_ancestor_term)
        for t in term_ids:
            term_id2clu[t] = (common_ancestor_term.get_id(), )
    return term_id2clu


def filter_clu_to_terms(term2clu):
    clu2term = invert_map(term2clu)
    for clu, terms in clu2term.items():
        if len(terms) == 1:
            del term2clu[terms.pop()]


def cover_with_onto_terms(model, onto, species_id2chebi_id, term_id2clu, ubiquitous_chebi_ids, r_ids_to_ignore=None):
    onto_updated = update_onto(onto, term_id2clu)
    if onto_updated:
        for clu, t_ids in invert_map(term_id2clu).items():
            if len(t_ids) == 1:
                del term_id2clu[t_ids.pop()]
            else:
                new_t_id2clu = cover_t_ids(model, species_id2chebi_id, ubiquitous_chebi_ids, t_ids, onto, clu,
                                           r_ids_to_ignore=r_ids_to_ignore)
                for t_id in t_ids:
                    if t_id in new_t_id2clu:
                        term_id2clu[t_id] = new_t_id2clu[t_id]
                    else:
                        del term_id2clu[t_id]
    return onto_updated


def maximization_step(model, onto, species_id2chebi_id, term_id2clu, ub_term_ids, unmapped_s_ids, r_ids_to_ignore=None):
    onto_updated = True
    while onto_updated:
        logging.info("  satisfying metabolite diversity...")
        term_id2clu = maximize(unmapped_s_ids, model, term_id2clu, species_id2chebi_id, ub_term_ids,
                               r_ids_to_ignore=r_ids_to_ignore)
        onto_updated = cover_with_onto_terms(model, onto, species_id2chebi_id, term_id2clu, ub_term_ids,
                                             r_ids_to_ignore=r_ids_to_ignore)


def find_term_clustering(model, chebi, species_id2chebi_id, unmapped_s_ids, ubiquitous_chebi_ids, r_ids_to_ignore=None):
    """
    Calculates a ChEBI term id clustering for the given model.
    :param model: libsbml.Model model of interest
    :param chebi: mod_sbml.onto.obo_ontology.Ontology ChEBI ontology
    :param species_id2chebi_id: dict {metabolite_id: ChEBI_term_id}
    :param unmapped_s_ids: set of ids of metabolite for which no ChEBI term was found
    :param ubiquitous_chebi_ids: set of ubiquitous ChEBI ids
    :param r_ids_to_ignore: (optional) ids of reactions whose stoichiometry preserving constraint can be ignores
    :return: dict {ChEBI_term_id: cluster}
    """
    if not ubiquitous_chebi_ids:
        ubiquitous_chebi_ids = set()
    chebi_ids = set(species_id2chebi_id.values()) - ubiquitous_chebi_ids

    logging.info("  aggressive metabolite grouping...")
    term_id2clu = cover_t_ids(model, species_id2chebi_id, ubiquitous_chebi_ids, chebi_ids, chebi,
                              r_ids_to_ignore=r_ids_to_ignore)
    chebi.trim({it[0] for it in term_id2clu.values()}, relationships=EQUIVALENT_RELATIONSHIPS)
    suggest_clusters(model, unmapped_s_ids, term_id2clu, species_id2chebi_id, ubiquitous_chebi_ids,
                     r_ids_to_ignore=r_ids_to_ignore)
    # filter_clu_to_terms(term_id2clu)
    # _log_clusters(term_id2clu, onto, model)

    maximization_step(model, chebi, species_id2chebi_id, term_id2clu, ubiquitous_chebi_ids, unmapped_s_ids,
                      r_ids_to_ignore=r_ids_to_ignore)
    # filter_clu_to_terms(term_id2clu)
    # _log_clusters(term_id2clu, onto, model)

    logging.info("  preserving stoichiometry...")
    fix_stoichiometry(model, term_id2clu, species_id2chebi_id, ubiquitous_chebi_ids, chebi,
                      r_ids_to_ignore=r_ids_to_ignore)
    # filter_clu_to_terms(term_id2clu)
    # _log_clusters(term_id2clu, onto, model)

    maximization_step(model, chebi, species_id2chebi_id, term_id2clu, ubiquitous_chebi_ids, unmapped_s_ids,
                      r_ids_to_ignore=r_ids_to_ignore)
    # filter_clu_to_terms(term_id2clu)
    # _log_clusters(term_id2clu, onto, model)

    return term_id2clu


def generalize_species(model, s_id2chebi_id, ub_s_ids, chebi, ub_chebi_ids, threshold=UBIQUITOUS_THRESHOLD,
                       r_ids_to_ignore=None):
    """
    Groups metabolites of the model into clusters.
    :param model: libsbml.Model model of interest
    :param s_id2chebi_id: dict {metabolite_id: ChEBI_term_id}
    :param ub_s_ids: collection of ubiquitous metabolite ids
    :param chebi: mod_sbml.onto.obo_ontology.Ontology ChEBI ontology
    :param ub_chebi_ids: collection of ubiquitous ChEBI term ids
    :param threshold: threshold for a metabolite to be considered as frequently participating in reactions
    and therefore ubiquitous
    :param r_ids_to_ignore: (optional) ids of reactions whose stoichiometry preserving constraint can be ignores
    :return:
    """
    unmapped_s_ids = {s.getId() for s in model.getListOfSpecies() if s.getId() not in s_id2chebi_id}
    term_id2clu = find_term_clustering(model, chebi, s_id2chebi_id, unmapped_s_ids, ub_chebi_ids,
                                       r_ids_to_ignore=r_ids_to_ignore)
    if term_id2clu:
        term_id2clu = select_representative_terms(term_id2clu, chebi)
        s_id2clu = compute_s_id2clu(unmapped_s_ids, model, s_id2chebi_id, term_id2clu)
        clu2s_ids = invert_map(s_id2clu)
        for s_ids in clu2s_ids.values():
            if len(s_ids) == 1:
                del s_id2clu[s_ids.pop()]
    else:
        s_id2clu = {}
    if not ub_s_ids:
        frequent_ch_ids = get_frequent_term_ids(model, threshold)
        ub_s_ids = select_metabolite_ids_by_term_ids(model, frequent_ch_ids) - set(s_id2clu.keys())
    # unmapped_s_ids = {s_id for s_id in unmapped_s_ids if s_id not in s_id2clu}
    # infer_clusters(model, unmapped_s_ids, s_id2clu, species_id2chebi_id, ub_chebi_ids)
    return s_id2clu, ub_s_ids


# =============================================================


def _log_clusters(term_id2clu, onto, model):
    clu2term = invert_map(term_id2clu)
    blueprint = []
    logging.info("-------------------\nquotient species sets:\n-------------------")
    for clu in sorted(clu2term.keys(), key=lambda k: -len(clu2term[k])):
        term_ids = clu2term[clu]
        if len(term_ids) == 1:
            continue
        blueprint.append(len(term_ids))
        logging.info("(%d)\t%s\n" % (len(term_ids), [
            onto.get_term(it).get_name() if onto.get_term(it) else model.getSpecies(
                it).getName() if model.getSpecies(it) else it for it in term_ids]))
    logging.info("Cluster sizes: %s\n-------------------\n\n" % sorted(blueprint, key=lambda s: -s))
