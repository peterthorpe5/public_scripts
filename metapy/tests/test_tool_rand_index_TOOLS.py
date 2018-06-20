#!/usr/bin/env python

"""Test program for a function in the module tools:
reformat_cdhit_clustrs
"""

import os
import shutil
from pycits.tools import NotExecutableError
from nose.tools import nottest, assert_equal
from pycits.Rand_index import pairwise_comparison_Rand
import sklearn


# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data", "Rand_index")
OUTDIR = os.path.join("tests", "test_out_Rand_index")
test_files = ['blastclust99', 'cdhit0.99',
              'perfect_map', 'swarm',
              'vsearch0.99']
FILE_LIST = []
for to_file in test_files:
    name = os.path.join("tests", "test_data", "Rand_index",
                        to_file)
    FILE_LIST.append(name)

RAND_RESULTS = []

# folder checking
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets",
                      "Rand_index",
                      "Rand_comparisons.txt")


def test_Rand_index_exec():
    """run the function Rand_index
    to compare the results against
    pregenerated results"""
    outfile = os.path.join(OUTDIR, "tests_Rand.txt")
    Rand_results = pairwise_comparison_Rand(sorted(FILE_LIST),
                                            outfile)

    with open(TARGET, "r") as target_fh:
        with open(outfile, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
