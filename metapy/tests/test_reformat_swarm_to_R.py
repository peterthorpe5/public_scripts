#!/usr/bin/env python

"""Test program for a function in the module tools:
reformat_swarm_cls
"""

import os
import shutil
from pycits.tools import NotExecutableError
from nose.tools import nottest, assert_equal
from pycits.tools import reformat_swarm_cls


# INPUT DATA LOCATION
# Note: We are reusing files for other tests here.
INDIR = os.path.join("tests", "test_targets", "swarm")
INFILE = "swarm.out"
OUTDIR = os.path.join("tests", "swarm_reform_for_R")
DBFILE = os.path.join("tests", "test_data",
                      "swarm",
                      "swarm_coded_with_abundance.fasta")

TEST_OUT = "swarm_result_Rout"

# folder checking
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets",
                      "swarm", "Format_for_R",
                      "swarm_tests_Rout")


def test_convert_exec():
    """run the function reformat_swarm_cls
    to compare the results against
    pregenerated results"""

    result_file = os.path.join(OUTDIR, TEST_OUT)

    # reformat_swarm_cls(swarm, db, db_and_reads, outfile)
    reformat_swarm_cls(os.path.join(INDIR, INFILE),
                       DBFILE,
                       DBFILE,
                       result_file,
                       True)

    with open(TARGET, "rt") as target_fh:
        data = target_fh.read().split("\n")
        with open(result_file, "r") as test_fh:
            result = test_fh.read().split("\n")

            assert_equal(result, data)
