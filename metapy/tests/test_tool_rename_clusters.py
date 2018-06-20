#!/usr/bin/env python

"""Test program for a function in the module tools: dereplicate_name
"""

import os
import shutil
import subprocess
from pycits.tools import NotExecutableError
import gzip
from nose.tools import nottest, assert_equal
from pycits.tools import parse_tab_file_get_clusters, coded_name_to_species


# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data", "rename_clusters")
OUTDIR = os.path.join("tests", "test_out_deduplicate_hash")
TEST_FILE = "swarm.out"
TEST_DB = "db_old_to_new_names.txt"
TEST_OUT = "swarm_test_nemaed.out"
OTU_DATABASE = os.path.join(INDIR, "ITS_database_NOT_confirmed" +
                            "_correct_last14bases_removed.fasta")


# folder checking
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "rename_clusters",
                      "swarm.outRENAMED")


def test_convert_exec():
    """Run parse_tab_file_get_clusters function from tool on test
    data and compare output to precomputed target"""
    # a cluster filen and  database of old name to new name already generated
    # need to compare renamed cluster file only.
    result_file = os.path.join(OUTDIR, TEST_OUT)
    db = os.path.join(INDIR, TEST_DB)
    # parse_tab_file_get_clusters(filename, database, out_file)
    parse_tab_file_get_clusters((os.path.join(INDIR, TEST_FILE)),
                                OTU_DATABASE, db, result_file,
                                True)
    # compare the precomputaed results versus these results.
    with open(TARGET, "rt") as target_fh:
        with open(result_file, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())

OTU_DATABASE = os.path.join(INDIR, "ITS_database_NOT_confirmed_correct" +
                            "_last14bases_removedabundance.fasta")
OTU_DATABASE = "./data/ITS_db_NOT_confirmed_for_swarm.fasta"
# run this again with another formatted db with abundance values.


def test_convert_exec():
    """Run parse_tab_file_get_clusters function from tool on test
    data and compare output to precomputed target - with abundance values"""
    # a cluster filen and  database of old name to new name already generated
    # need to compare renamed cluster file only.
    result_file = os.path.join(OUTDIR, TEST_OUT)
    db = os.path.join(INDIR, TEST_DB)
    # parse_tab_file_get_clusters(filename, database, out_file)
    parse_tab_file_get_clusters((os.path.join(INDIR, TEST_FILE)),
                                OTU_DATABASE, db, result_file,
                                True)
    # compare the precomputaed results versus these results.
    with open(TARGET, "rt") as target_fh:
        with open(result_file, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
