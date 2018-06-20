#!/usr/bin/env python

"""Test program for a function in the module tools: filter_sam_file
"""

import os
import shutil
import subprocess
from pycits.tools import NotExecutableError
import pysam
import gzip
from nose.tools import nottest, assert_equal
from pycits.tools import filter_sam_file


# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data", "Samfile")
OUTDIR = os.path.join("tests", "test_out_Samfile")
TEST_FILE = "tests_filter.sam"
TEST_OUT = "generated.sam"


# folder checking
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# TARGET OUTPUT DATA
TARGET_SAMFILE = os.path.join("tests", "test_targets", "Samfile",
                              "sam.target")
TARGET_CIGLIST = os.path.join("tests", "test_targets", "Samfile",
                              "cig_list.txt")

CIGLIST = [[(0, 171)], [(0, 171)], [(0, 194)],
           [(0, 213)], [(0, 213)], [(0, 218)]]


def test_convert_exec():
    """run the function to compare the results against
    pregenerated results"""

    result_file = os.path.join(OUTDIR, TEST_OUT)
    # parse_tab_file_get_clusters(filename, database, out_file)
    cig_list, matches = filter_sam_file(os.path.join(INDIR, TEST_FILE),
                                        result_file)
    data = []
    for entry in matches:
        data.append(entry.rstrip())

    assert_equal(cig_list, CIGLIST)
    if not os.path.isfile(TARGET_SAMFILE):
        print ("cant find the target")
    premade = []
    with open(TARGET_SAMFILE, "rt") as target_fh:
        for line in target_fh:
            premade.append(line.rstrip())

    assert_equal(premade, data)
