#!/usr/bin/env python
# author: Peter Thorpe October 2019.
# U of St Andrew
# script to test the from HGT_modules.parsenodes import parse_NCBI_nodes_tab_file
# pip install nose!!!

import os
from nose.tools import nottest, assert_equal, assert_not_equal
import unittest
from HGT_modules.parseblast import parse_blast_line


# INPUT DATA LOCATION
INPUT = os.path.join("tests", "inputs")

BLAST_LINE = "GPALN_014324-T1 gi|612340452|gb|AHW98762.1|     93.5    336     21      1       1       335     1       336     2.2e-186        661.0beta-1,4-endoglucanase precursor [Globodera rostochiensis]       31243   Globodera rostochiensis         Eukaryota"

def test_parse_NCBI_nodes_tab_file():
    "test parse_blast_line make sure taxid = tax id "
    query_name, percentage_identity, Evalue, bit_score, \
           description, tax_id, species_sci, species_common,\
           kingdom = parse_blast_line(BLAST_LINE, 14)
    assert_equal(tax_id, "31243")

def test_parse_NCBI_nodes_tab_file():
    "test parse_blast_line test name is not tax id "
    query_name, percentage_identity, Evalue, bit_score, \
           description, tax_id, species_sci, species_common,\
           kingdom = parse_blast_line(BLAST_LINE, 14)
    assertNotEqual(query_name, "31243")

def test_parse_NCBI_nodes_tab_file():
    "test parse_blast_line diff tax_id_column "
    query_name, percentage_identity, Evalue, bit_score, \
           description, tax_id, species_sci, species_common,\
           kingdom = parse_blast_line(BLAST_LINE, 13)
    assert_not_equal(tax_id, "31243")
