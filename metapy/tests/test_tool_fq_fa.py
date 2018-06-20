#!/usr/bin/env python

"""Test program for a function in the module tools: convert_fq_to_fa
"""

import gzip
import os
import shutil
import subprocess

from nose.tools import nottest, assert_equal
from pycits.tools import convert_fq_to_fa
from pycits.tools import NotExecutableError


# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out", "tools")
INFILE = os.path.join(INDIR, "pear_test.assembled.fastq.gz")
OUTFILE = os.path.join(OUTDIR, "converted_fq_fa.fasta")

# folder checking
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "tools",
                      "pear_tests_FASTA_convert.fasta.gz")


# Create output directory tree
def setup():
    """Set up test fixtures"""
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)


def test_convert_exec():
    """Run convert_fq_to_fa.py on test data and compare output to
    precomputed target
    """
    convert_fq_to_fa(INFILE, OUTFILE)

    with gzip.open(TARGET, "rt") as target_fh:
        with open(OUTFILE, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
