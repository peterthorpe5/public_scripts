#!/usr/bin/env python

"""Test program for a function in the module tools:
check_OTU_db_abundance_val
"""

import os
import shutil
import subprocess
from pycits.tools import NotExecutableError
import gzip
from nose.tools import nottest, assert_equal
from pycits.tools import check_OTU_db_abundance_val
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data")
TEST_FILE = "test_db_for_abundance_checking.fasta"
TEST_OUT = "test_db_for_abundance_checkingabundance.fasta"

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "db_abundance_checking",
                      "target_db_for_abundance_checking.fasta")


def test_convert_exec():
    """Run check_OTU_db_abundance_val on test data and compare
    output to
    precomputed target
    """
    # a fasta file and a database of old name to new name is generated
    # need to compare both of these.
    result_file = os.path.join(INDIR, TEST_OUT)
    check_OTU_db_abundance_val((os.path.join(INDIR, TEST_FILE)))

    with open(TARGET, "r") as target_fh:
        with open(result_file, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
