#!/usr/bin/env python
# encoding: utf-8

import getopt
import logging
from os import listdir
import sys

import libsbml
import openpyxl
import openpyxl.styles
from openpyxl.styles import Style, Font
from openpyxl.styles.colors import Color

from mod_sbml.annotation.chebi.chebi_annotator import get_term
from mod_sbml.sbml.sbml_manager import get_gene_association
from sbml_generalization.merge.model_merger import merge_models
from sbml_generalization.generalization.sbml_generalizer import generalize_model
from mod_sbml.annotation.chebi.chebi_serializer import get_chebi
from mod_sbml.onto import parse_simple
from mod_sbml.utils.misc import invert_map
from mod_sbml.serialization import get_sbml_r_formula

__author__ = 'anna'

##
# runner module generalizes the model.
# usage: main.py --model model.xml --chebi chebi.obo --verbose
##

help_message = '''
Generalizes the model.
usage: main.py --model model.xml --verbose
'''

BASIC_STYLE = Style(font=Font(color=openpyxl.styles.colors.BLACK, sz=8))
BOLD_STYLE = Style(font=Font(bold=True, sz=8))
HEADER_STYLE = BOLD_STYLE


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        chebi, in_sbml, in_path, out_sbml, groups_sbml, merged_sbml, verbose, log_file = process_args(argv)
        if verbose:
            logging.basicConfig(level=logging.INFO)
        logging.info("parsing ChEBI...")
        ontology = parse_simple(chebi)
        if not in_sbml and in_path:
            in_sbml_list = ['%s/%s' % (in_path, f) for f in listdir(in_path)
                            if f.find(".xml") != -1 or f.find(".sbml") != -1]
            for sbml in in_sbml_list:
                doc = libsbml.SBMLReader().readSBML(sbml)
                model = doc.getModel()
                print(sbml, model.getNumSpecies(), model.getNumReactions())
            merge_models(in_sbml_list, merged_sbml)
            in_sbml = merged_sbml
        r_id2clu, s_id2clu, _, _ = generalize_model(groups_sbml, out_sbml, in_sbml, ontology, ub_chebi_ids={'chebi:ch'})
        if in_path:
            serialize_generalization(r_id2clu, s_id2clu, groups_sbml, ontology, '%s/generalized.xlsx' % in_path)
    except Usage, err:
        logging.error(sys.argv[0].split("/")[-1] + ": " + str(err.msg))
        logging.error(sys.stderr, "\t for help use --help")
        return 2


def add_values(ws, row, col, values, style=BASIC_STYLE):
    j = col
    for v in values:
        ws.cell(row=row, column=j).value = v
        ws.cell(row=row, column=j).style = style
        j += 1


def serialize_generalization(r_id2clu, s_id2clu, sbml, chebi, path):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    groups_plugin = model.getPlugin("groups")
    clu2r_ids, clu2s_ids = invert_map(r_id2clu), invert_map(s_id2clu)
    wb = openpyxl.Workbook()
    ws = wb.create_sheet(0, "Metabolite Groups")
    row = 1
    add_values(ws, row, 1, ["Group Id", "Group Name", "Group CHEBI", "Id", "Name", "Compartment", "CHEBI"], HEADER_STYLE)
    row += 1
    processed_s_ids = set()
    for (g_id, ch_term), s_ids in sorted(clu2s_ids.iteritems(), key=lambda ((g_id, _), s_ids): g_id):
        group = groups_plugin.getGroup(g_id)
        add_values(ws, row, 1, [g_id, group.getName(), ch_term.get_id() if ch_term else ''])
        for s_id in sorted(s_ids, key=lambda s_id: s_id[s_id.find('__'):]):
            species = model.getSpecies(s_id)
            ch_term = get_term(species, chebi)
            add_values(ws, row, 4, [s_id, species.getName(), model.getCompartment(species.getCompartment()).getName(),
                                    ch_term.get_id() if ch_term else ''])
            row += 1
        processed_s_ids |= s_ids
    ws = wb.create_sheet(1, "Ungrouped metabolites")
    row = 1
    add_values(ws, row, 1, ["Id", "Name", "Compartment", "CHEBI"], HEADER_STYLE)
    row += 1
    unm_l = 0
    for species in sorted(model.getListOfSpecies(), key=lambda s: s.getId()[s.getId().find('__'):]):
        if species.getId() not in processed_s_ids:
            ch_term = get_term(species, chebi)
            add_values(ws, row, 1, [species.getId(), species.getName(),
                                    model.getCompartment(species.getCompartment()).getName(),
                                    ch_term.get_id() if ch_term else ''])
            row += 1
            unm_l += 1
    print unm_l

    ws = wb.create_sheet(2, "Reaction Groups")
    row = 1
    add_values(ws, row, 1, ["Group Id", "Id", "Name", "Formula", "Gene Association"], HEADER_STYLE)
    row += 1
    processed_r_ids = set()
    for (g_id, g_name), r_ids in sorted(clu2r_ids.iteritems(), key=lambda ((g_id, g_name), _): g_id):
        add_values(ws, row, 1, [g_id])
        for r_id in sorted(r_ids, key=lambda r_id: r_id[r_id.find('__'):]):
            r = model.getReaction(r_id)
            add_values(ws, row, 2, [r_id, r.getName(), get_sbml_r_formula(model, r, show_compartments=True, show_metabolite_ids=True),
                                    get_gene_association(r)])
            row += 1
        processed_r_ids |= r_ids
    ws = wb.create_sheet(3, "Ungrouped reactions")
    row = 1
    add_values(ws, row, 1, ["Id", "Name", "Formula", "Gene Association"], HEADER_STYLE)
    row += 1
    unm_l = 0
    for r in sorted(model.getListOfReactions(), key=lambda r: r.getId()[r.getId().find('__'):]):
        if r.getId() not in processed_r_ids:
            add_values(ws, row, 1, [r.getId(), r.getName(), get_sbml_r_formula(model, r, show_compartments=True, show_metabolite_ids=True),
                                    get_gene_association(r)])
            row += 1
            unm_l += 1
    print unm_l

    wb.save(path)


def generate_out_sbml_name(in_sbml, in_path, out_sbml, groups_sbml, merged_sbml):
    if not in_sbml and not in_path:
        raise Usage(help_message)
    if in_sbml:
        extension = in_sbml.find(".xml")
        if extension == -1:
            extension = in_sbml.find(".sbml")
        if extension == -1:
            extension = len(in_sbml)
        prefix = in_sbml[:extension]
    else:
        prefix = '{0}/merged'.format(in_path)
        if not merged_sbml:
            merged_sbml = "{0}.xml".format(prefix)
    if not out_sbml:
        out_sbml = "{0}_generalized.xml".format(prefix)
    if not groups_sbml:
        groups_sbml = "{0}_with_groups.xml".format(prefix)
    return out_sbml, groups_sbml, merged_sbml


def process_args(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "m:c:h:o:v:l:p:g:b",
                                   ["help", "model=", "chebi=", "outmodel=", "verbose", "log=", "inpath=", "groupsmodel=",
                                    "mergedmodel="])
    except getopt.error, msg:
        raise Usage(msg)
    in_path, in_sbml, groups_sbml, merged_sbml, chebi, out_sbml, verbose, log_file = \
        None, None, None, None, None, None, False, None
    # option processing
    for option, value in opts:
        if option in ("-h", "--help"):
            raise Usage(help_message)
        if option in ("-m", "--model"):
            in_sbml = value
        if option in ("-c", "--chebi"):
            chebi = value
        if option in ("-o", "--outmodel"):
            out_sbml = value
        if option in ("-g", "--groupsmodel"):
            groups_sbml = value
        if option in ("-b", "--mergedmodel"):
            merged_sbml = value
        if option in ("-v", "--verbose"):
            verbose = True
        if option in ("-l", "--log"):
            log_file = value
        if option in ("-p", "--inpath"):
            in_path = value
    out_sbml, groups_sbml, merged_sbml = generate_out_sbml_name(in_sbml, in_path, out_sbml, groups_sbml, merged_sbml)
    if not chebi:
        chebi = get_chebi()
    if not in_sbml and not in_path:
        raise Usage(help_message)
    return chebi, in_sbml, in_path, out_sbml, groups_sbml, merged_sbml, verbose, log_file


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

if __name__ == "__main__":
    sys.exit(main())