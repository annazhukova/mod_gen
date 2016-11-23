from collections import defaultdict
from itertools import chain
import threading

from sbml_generalization.generalization.vertical_key import is_reactant

__author__ = 'anna'

max_lock = threading.RLock()


def merge_based_on_neighbours(lst):
    new_lst = []
    for neighbours, terms in lst:
        neighbours = set(neighbours)
        to_remove = []
        for (new_neighbours, new_terms) in new_lst:
            if neighbours & new_neighbours:
                neighbours |= new_neighbours
                terms |= new_terms
                to_remove.append((new_neighbours, new_terms))
        new_lst = [it for it in new_lst if not it in to_remove] + [(neighbours, terms)]
    return new_lst


class MaximizingThread(threading.Thread):
    def __init__(self, model, term_ids, species_id2term_id, clu, term_id2clu, s_id2clu,
                 ubiquitous_chebi_ids, r_id2clu, r_ids_to_ignore=None):
        threading.Thread.__init__(self)
        self.model = model
        self.term_ids = term_ids
        self.species_id2term_id = species_id2term_id
        self.clu = clu
        self.term_id2clu = term_id2clu
        self.s_id2clu = s_id2clu
        self.ubiquitous_chebi_ids = ubiquitous_chebi_ids
        self.r_id2clu = r_id2clu
        self.r_ids_to_ignore = r_ids_to_ignore

    def run(self):
        neighbours2term_ids = defaultdict(set)
        neighbourless_terms = set()
        t_id2rs = defaultdict(list)
        for r in (r for r in self.model.getListOfReactions() if r.getNumReactants() + r.getNumProducts() > 2):
            if self.r_ids_to_ignore and r.getId() in self.r_ids_to_ignore:
                continue
            for s_id in chain((species_ref.getSpecies() for species_ref in r.getListOfReactants()),
                              (species_ref.getSpecies() for species_ref in r.getListOfProducts())):
                if s_id in self.species_id2term_id:
                    t_id2rs[self.species_id2term_id[s_id]].append(r)
                else:
                    t_id2rs[s_id].append(r)
        for t_id in self.term_ids:
            neighbours = {
                ("in"
                 if is_reactant(self.model, t_id, r, self.s_id2clu, self.species_id2term_id, self.ubiquitous_chebi_ids)
                 else "out",
                 self.r_id2clu[r.getId()]) for r in t_id2rs[t_id]}
            if neighbours:
                key = tuple(sorted(neighbours))
                neighbours2term_ids[key].add(t_id)
            else:
                neighbourless_terms.add(t_id)
        new_lst = merge_based_on_neighbours(neighbours2term_ids.items())
        i = 0
        if len(new_lst) > 1:
            for neighbours, term_ids in new_lst:
                n_clu = self.clu + (i,)
                i += 1
                with max_lock:
                    for t in term_ids:
                        self.term_id2clu[t] = n_clu
        with max_lock:
            for t in neighbourless_terms:
                self.term_id2clu[t] = self.clu + (i,)
                i += 1
