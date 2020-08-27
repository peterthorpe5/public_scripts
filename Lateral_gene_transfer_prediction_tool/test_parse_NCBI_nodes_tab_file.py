#!/usr/bin/env python
# author: Peter Thorpe October 2019.
# U of St Andrew
# script to test the from HGT_modules.parsenodes import parse_NCBI_nodes_tab_file
# pip install nose!!!

import os
from nose.tools import nottest, assert_equal
import unittest
from HGT_modules.parsenodes import parse_NCBI_nodes_tab_file


# INPUT DATA LOCATION
INPUT = os.path.join("tests", "inputs")

print("here")
test_dict = parse_NCBI_nodes_tab_file(INPUT)
print(test_dict)


def test_parse_NCBI_nodes_tab_file():
    "test parse_NCBI_nodes_tab_file 1 "
    test_dict = parse_NCBI_nodes_tab_file(os.path.join("tests",
                                                       "inputs"))
    print(test_dict)
##    assert_equal(intergenic_region,
##                 "GTGAG")
##
##def test_intergenic_region4():
##    "test intergenic_region 4 "
##    Genome_seq_record = Genome_sequence["pathogens_Gpal_scaffold_1"]
##    intergenic_region = get_len_upstream(Genome_seq_record.seq,
##                                         4, "+")
##    assert_equal(intergenic_region,
##                 "TGAG")
##
##def test_intergenic_region3():
##    "test intergenic_region 3 "
##    Genome_seq_record = Genome_sequence["pathogens_Gpal_scaffold_1"]
##    intergenic_region = get_len_upstream(Genome_seq_record.seq,
##                                         3, "+")
##    assert_equal(intergenic_region,
##                 "GAG")
##    
#################################################################
### - ve coded test
##
##def test_negative_region3():
##    "test Rev intergenic_region 3 "
##    Genome_seq_record = Genome_sequence["pathogens_Gpal_scaffold_2"]
##    intergenic_region = get_len_upstream(Genome_seq_record.seq,
##                                         3, "-")
##    # this has already been reverse complemented, so we want the
##    # seq at the "end" of the ROI. 
##    assert_equal(intergenic_region,
##                 "NNN")
##
##
##def test_negative_region3_short_seq():
##    "test Rev intergenic_9a "
##    Genome_seq_record = Genome_sequence["pathogens_Gpal_scaffold_3"]
##    intergenic_region = get_len_upstream(Genome_seq_record.seq,
##                                         8, "-")
##    # this has already been reverse complemented, so we want the
##    # seq at the "end" of the ROI.
##    assert_equal(intergenic_region,
##                 "ATGNNATG")
##
##
##def test_positive_region3_short_seq():
##    "test Rev  9b "
##    Genome_seq_record = Genome_sequence["pathogens_Gpal_scaffold_3"]
##    intergenic_region = get_len_upstream(Genome_seq_record.seq,
##                                         8, "+")
##    # this has already been reverse complemented, so we want the
##    # seq at the "end" of the ROI. 
##    assert_equal(intergenic_region,
##                 "ATGNNAT")
##
##
##def test_positive_region3_short_seq():
##    "test Rev  9c "
##    Genome_seq_record = Genome_sequence["pathogens_Gpal_scaffold_3"]
##    intergenic_region = get_len_upstream(Genome_seq_record.seq,
##                                         15, "+")
## 
##    assert_equal(intergenic_region,
##                 "ATGNNATG")
##
##
##def test_positive_region3_short_seq():
##    "test Rev  9d "
##    Genome_seq_record = Genome_sequence["pathogens_Gpal_scaffold_3"]
##    intergenic_region = get_len_upstream(Genome_seq_record.seq,
##                                         14, "+") 
##    assert_equal(intergenic_region,
##                 "ATGNNATG") 
