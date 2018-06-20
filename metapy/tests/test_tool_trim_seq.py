#!/usr/bin/env python

"""Test program for a function in the module tools: trim_seq
"""

import os
import shutil
import subprocess
from pycits.tools import NotExecutableError
import gzip
from nose.tools import nottest, assert_equal
from pycits.tools import trim_seq


# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_trim_seq")
TEST_FILE = "test_for_trim_function.fasta"
TEST_OUT = "test_for_trim_function_out_L53_R0.fasta"
# folder checking
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "trim_seq",
                      "test_for_trim_function_out_L53_R0.fasta")


def test_convert_exec():
    """Run trim_seq on test data and compare output to
    precomputed target left 53 right 0 - default values
    """
    result_file = os.path.join(OUTDIR, TEST_OUT)
    trim_seq((os.path.join(INDIR, TEST_FILE)), result_file)

    with open(TARGET, "rt") as target_fh:
        with open(result_file, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())


# test again with different Left and Right value
TEST_OUT_L15 = "test_for_trim_function_out_L15_R10.fasta"
# TARGET OUTPUT DATA - NEW
TARGET_L15 = os.path.join("tests", "test_targets", "trim_seq",
                          "test_for_trim_function_out_L15_R10.fasta")


def test_convert_exec():
    """Run convert_trim_seq on test data and compare output to
    precomputed target left 15 right 10
    """
    result_file = os.path.join(OUTDIR, TEST_OUT_L15)
    trim_seq((os.path.join(INDIR, TEST_FILE)), result_file, 15, 10)

    with open(TARGET_L15, "rt") as target_fh:
        with open(result_file, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
