#!/usr/bin/env python

"""Test program for a function in the module tools: dereplicate_name
"""

import os
import shutil
import subprocess
from pycits.tools import NotExecutableError
import gzip
from nose.tools import nottest, assert_equal
from pycits.tools import dereplicate_name, deduplicate_and_rename
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_deduplicate_hash")
TEST_FILE = "test_for_dedup_rename_hash.fasta"
TEST_OUT = "test_OUTPUT_dedup_rename_hash.fasta"
TEST_DB_OUT = "test_DB.txt"
# folder checking
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "deduplicate_rename",
                      "tests_deduplicated_hash.fasta")
TARGET_DB = os.path.join("tests", "test_targets", "deduplicate_rename",
                         "database_tests.txt")


def return_id_seq(fasta, my_list):
    """function to populate a list for the purpose of comparison.
    you cannot directly compare a bioseqrecord.
    """
    for seq_record in SeqIO.parse(fasta, 'fasta'):
        data = "%s%s" % (seq_record.id, seq_record.seq)
        my_list.append(data)
    return sorted(my_list)


def populate_my_db_list(db_file, db_list):
    """function to put the contents of the db file in a list,
    so it can be sorted and compared"""
    with open(db_file, "r") as fh:
        for line in fh:
            db_list.append(line)
    return sorted(db_list)


def test_convert_exec():
    """Run dereplicate_name on test data and compare output to
    precomputed target l
    """
    # a fasta file and a database of old name to new name is generated
    # need to compare both of these.
    result_file = os.path.join(OUTDIR, TEST_OUT)
    result_db = os.path.join(OUTDIR, TEST_DB_OUT)
    dereplicate_name((os.path.join(INDIR, TEST_FILE)), result_db, result_file)

    # list to populate with the data for sorting and comparing
    tests_fa_list = list()
    result_db_list = list()
    target_fa_list = list()
    target_db_list = list()
    # call function to pouplate the relavant lists
    result_db_list = populate_my_db_list(result_db, result_db_list)
    target_db_list = populate_my_db_list(TARGET_DB, target_db_list)

    tests_fa_list = return_id_seq(result_file, tests_fa_list)
    target_fa_list = return_id_seq(TARGET, target_fa_list)
    # these dont come back in the same order so, need to sort them
    assert_equal(tests_fa_list, target_fa_list)

    assert_equal(result_db_list, target_db_list)
