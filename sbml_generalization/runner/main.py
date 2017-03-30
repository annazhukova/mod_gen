#!/usr/bin/env python
# encoding: utf-8

import logging
import os

from mod_sbml.annotation.chebi.chebi_serializer import get_chebi
from mod_sbml.onto import parse_simple
from sbml_generalization.generalization.sbml_generalizer import generalize_model

__author__ = 'anna'


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generalizes an SBML model.")
    parser.add_argument('--model', required=True, type=str,
                        help="input model in SBML format")
    parser.add_argument('--output_model', default=None, type=str,
                        help="path to the output generalized model in SBML format")
    parser.add_argument('--groups_model', default=None, type=str,
                        help="path to the output model in SBML format with groups extension to encode similar elements")
    parser.add_argument('--verbose', action="store_true", help="print logging information")
    parser.add_argument('--log', default=None, help="a log file")
    params = parser.parse_args()

    prefix = os.path.splitext(params.model)[0]
    if not params.output_model:
        params.output_model = "%s_generalized.xml" % prefix
    if not params.groups_model:
        params.groups_model = "%s_with_groups.xml" % prefix

    if params.verbose:
        logging.basicConfig(level=logging.INFO)

    logging.info("parsing ChEBI...")
    ontology = parse_simple(get_chebi())
    r_id2clu, s_id2clu, _, _ = generalize_model(params.model, ontology, params.groups_model, params.output_model,
                                                ub_chebi_ids={'chebi:ch'})
