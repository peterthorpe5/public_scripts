#!/usr/bin/env python

"""Test program for a function in the module tools:
reformat_sam_clusters
"""

import os
import shutil
from pycits.tools import NotExecutableError
from nose.tools import nottest, assert_equal
from pycits.tools import reformat_sam_clusters


# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data", "bowtie",
                     "for_R_format")
OUTDIR = os.path.join("tests", "sam_reform_for_R")
TEST_FILE_DB = "tests_db.fasta"
TEST_FILE_SAM = "tests_filter.sam"

TEST_OUT = "bowtie_result_Rout"

# folder checking
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets",
                      "Samfile", "reformat_for_R",
                      "bowtie_tests_Rout")


def test_convert_exec():
    """run the function reformat_sam_clusters
    to compare the results against
    pregenerated results"""

    result_file = os.path.join(OUTDIR, TEST_OUT)

    # reformat_sam_clusters(sam, db_and_reads, outfile)
    reformat_sam_clusters(os.path.join(INDIR, TEST_FILE_SAM),
                          os.path.join(INDIR, TEST_FILE_DB),
                          result_file)

    with open(TARGET, "rt") as target_fh:
        data = target_fh.read().split("\n")
        sorted_data = sorted(data)

        with open(result_file, "r") as test_fh:
            result = test_fh.read().split("\n")
            sorted_result = sorted(result)

            assert_equal(sorted_result, sorted_data)

test_convert_exec()
