#!/usr/bin/env python

"""Test program for a function in the module tools:
reformat_cdhit_clustrs
"""

import os
import shutil
from pycits.tools import NotExecutableError
from nose.tools import nottest, assert_equal
from pycits.tools import reformat_cdhit_clustrs


# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data", "reformat_cdhit")
OUTDIR = os.path.join("tests", "cdhit_reformat")
TEST_FILE = "test.clstr"
TEST_OUT_1_LINE_PER = "cdhit_1_line_per.clstr"
TEST_OUT_R_FORMAT = "cdhit_r_format.clstr"

# folder checking
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# TARGET OUTPUT DATA
TARGET_1_LINE_PER = os.path.join("tests", "test_targets",
                                 "reformat_cdhit",
                                 "test_one_line")
TARGET_R_FORMAT = os.path.join("tests", "test_targets",
                               "reformat_cdhit",
                               "tests_Rout")


def test_convert_exec():
    """run the function reformat_cdhit_clustrs
    to compare the results against
    pregenerated results"""

    result_file_1_line_per = os.path.join(OUTDIR, TEST_OUT_1_LINE_PER)
    result_file_R = os.path.join(OUTDIR, TEST_OUT_R_FORMAT)

    # reformat_cdhit_clustrs(clustr, outfile, out_R)
    reformat_cdhit_clustrs(os.path.join(INDIR, TEST_FILE),
                           result_file_1_line_per, result_file_R)

    with open(TARGET_1_LINE_PER, "rt") as target_fh:
        with open(result_file_1_line_per, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())

    with open(TARGET_R_FORMAT, "rt") as target_fh:
        with open(result_file_R, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
