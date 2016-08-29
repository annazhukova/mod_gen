import logging

import libsbml
from mod_sbml.annotation.chebi.chebi_annotator import annotate_metabolites

from mod_sbml.annotation.gene_ontology.go_annotator import get_go_id, annotate_compartments
from mod_sbml.annotation.gene_ontology.go_serializer import get_go
from mod_sbml.sbml.compartment.compartment_manager import need_boundary_compartment, \
    separate_boundary_metabolites
from mod_sbml.annotation.rdf_annotation_helper import get_qualifier_values, add_annotation
from mod_sbml.onto import parse_simple
from mod_sbml.annotation.chebi.chebi_serializer import get_chebi
from sbml_generalization.sbml.sbml_helper import set_consistency_level


__author__ = 'anna'


CYTOPLASM = 'go:0005737'
CYTOSOL = 'go:0005829'


def update_model_element_ids(m_id, model, go2c_id, go, chebi):
    id2id = {}
    if need_boundary_compartment(model):
        separate_boundary_metabolites(model)
    annotate_metabolites(model, chebi)
    annotate_compartments(model, go)

    for c in model.getListOfCompartments():
        c_id = c.getId()
        go_id = get_go_id(c)
        if not go_id:
            go_id = c.getName()
        if go_id:
            go_id = go_id.lower()
            if go_id == CYTOSOL and (CYTOSOL not in go2c_id) and (CYTOPLASM in go2c_id):
                go_id = CYTOPLASM
            elif go_id == CYTOPLASM and (CYTOSOL in go2c_id) and (CYTOPLASM not in go2c_id):
                go_id = CYTOSOL
            if go_id not in go2c_id:
                go2c_id[go_id] = "c_%s__%s" % (m_id, c_id)
            new_c_id = go2c_id[go_id]
        else:
            new_c_id = "c_%s__%s" % (m_id, c_id)
        c.setId(new_c_id)
        id2id[c_id] = new_c_id

    for c in model.getListOfCompartments():
        if c.getOutside():
            c.setOutside(id2id[c.getOutside()])

    for s in model.getListOfSpecies():
        if s.getCompartment():
            s.setCompartment(id2id[s.getCompartment()])
        old_id = s.getId()
        new_id = "s_%s__%s" % (m_id, old_id)
        s.setId(new_id)
        id2id[old_id] = new_id
        if s.getSpeciesType():
            s.unsetSpeciesType()

    for r in model.getListOfReactions():
        old_id = r.getId()
        new_id = "r_%s__%s" % (m_id, old_id)
        r.setId(new_id)
        id2id[old_id] = new_id
        if r.getCompartment():
            r.setCompartment(id2id[r.getCompartment()])
        for s_ref in r.getListOfReactants():
            s_ref.setSpecies(id2id[s_ref.getSpecies()])
        for s_ref in r.getListOfProducts():
            s_ref.setSpecies(id2id[s_ref.getSpecies()])
        for s_ref in r.getListOfModifiers():
            s_ref.setSpecies(id2id[s_ref.getSpecies()])


def get_model_id(i, m_ids, model):
    m_id = ''.join(e for e in model.getId() if e.isalnum()) if model.getId() else "m"
    if m_id in m_ids:
        while "m_%s_%d" % (m_id, i) in m_ids:
            i += 1
        m_id = "m_%s_%d" % (m_id, i)
        model.setId(m_id)
    m_ids.add(m_id)
    return m_id


def copy_reaction(e, model):
    new_e = model.createReaction()
    new_e.setId(e.getId())
    new_e.setName(e.getName())
    new_e.setCompartment(e.getCompartment())
    new_e.setReversible(e.getReversible())
    new_e.setKineticLaw(e.getKineticLaw())
    for s in e.getListOfReactants():
        sr = new_e.createReactant()
        sr.setSpecies(s.getSpecies())
        sr.setStoichiometry(s.getStoichiometry())
    for s in e.getListOfProducts():
        sr = new_e.createProduct()
        sr.setSpecies(s.getSpecies())
        sr.setStoichiometry(s.getStoichiometry())
    for s in e.getListOfModifiers():
        sr = new_e.createModifier()
        sr.setSpecies(s.getSpecies())
        sr.setStoichiometry(s.getStoichiometry())
    new_e.setNotes(e.getNotes())
    return new_e


def copy_species(e, model):
    new_e = model.createSpecies()
    new_e.setId(e.getId())
    new_e.setName(e.getName())
    new_e.setCompartment(e.getCompartment())
    new_e.setBoundaryCondition(e.getBoundaryCondition())
    new_e.setNotes(e.getNotes())
    return new_e


def copy_compartment(e, model):
    new_e = model.createCompartment()
    new_e.setId(e.getId())
    new_e.setName(e.getName())
    new_e.setOutside(e.getOutside())
    new_e.setNotes(e.getNotes())
    return new_e


def merge_models(in_sbml_list, out_sbml):
    if not in_sbml_list:
        raise ValueError('Provide SBML models to be merged')
    go = parse_simple(get_go())
    chebi = parse_simple(get_chebi())
    i = 0
    model_ids = set()
    go2c_id = {}

    doc = libsbml.SBMLDocument(2, 4)
    model = doc.createModel()
    model.setId('m_merged')
    m_c_ids = set()

    for o_sbml in in_sbml_list:
        o_doc = libsbml.SBMLReader().readSBML(o_sbml)
        set_consistency_level(o_doc)
        o_doc.checkL2v4Compatibility()
        o_doc.setLevelAndVersion(2, 4, False, True)
        o_model = o_doc.getModel()
        logging.info("Processing %s" % o_sbml)
        model_id = get_model_id(i, model_ids, o_model)

        update_model_element_ids(model_id, o_model, go2c_id, go, chebi)
        for e in o_model.getListOfCompartments():
            c_id = e.getId()
            if c_id not in m_c_ids:
                if model.addCompartment(e):
                    copy_compartment(e, model)
                m_c_ids.add(c_id)
        for e in o_model.getListOfSpecies():
            if model.getSpecies(e.getId()):
                continue
            if model.addSpecies(e):
                copy_species(e, model)
        for e in o_model.getListOfReactions():
            if model.addReaction(e):
                copy_reaction(e, model)

    libsbml.writeSBMLToFile(doc, out_sbml)