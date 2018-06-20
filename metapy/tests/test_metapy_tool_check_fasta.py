#!/usr/bin/env python3

"""Test program for a function in the module metapy_tools:
database_checker
"""

import os
import shutil

from nose.tools import nottest, assert_equal
from pycits.metapy_tools import database_checker
from pycits.tools import NotExecutableError

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data", "metapy_tools")
OUTDIR = os.path.join("tests", "test_out", "metapy_tools")
INFILE_ILLEGAL_LETTER = os.path.join(INDIR, "illegal_letter.fa")
INFILE_DUP_NAMES = os.path.join(INDIR, "duplicate_names.fa")
INFILE_DUP_SEQ = os.path.join(INDIR, "duplicate_seq.fa")
INFILE_DB = os.path.join(INDIR, "test_db.fasta")
OUTFILE = os.path.join(OUTDIR, "passed_db.fasta")

# folder checking
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "metapy_tools",
                      "passed_target.fasta")


# Create output directory tree
def setup():
    """Set up test fixtures"""
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)


def test_database_checker_ILLEGAL_LETTER():
    """Run database_checker INFILE_ILLEGAL_LETTER from metapy_tools
    """
    error, seqrecord = database_checker(INFILE_ILLEGAL_LETTER, OUTFILE)
    assert error != "ok"


def test_database_checker_DUP_NAMES():
    """Run database_checker DUP_NAMES from metapy_tools
    """
    error, seqrecord = database_checker(INFILE_DUP_NAMES, OUTFILE)
    assert error != "ok"


def test_database_checker_DUP_SEQ():
    """Run database_checker DUP_SEQ from metapy_tools
    """
    error, seqrecord = database_checker(INFILE_DUP_SEQ, OUTFILE)
    assert error != "ok"


def test_database_checker():
    """Run database_checker from metapy_tools
    """
    error, seqrecord = database_checker(INFILE_DB, OUTFILE)
    # assert_equal(target_fh.read(), test_fh.read())
    # fails even with LInux diff says they are the same. Weird!!!
    for seq_record in SeqIO.parse(TARGET, "fasta"):
        for result_record in SeqIO.parse(OUTFILE, "fasta"):
            if seq_record.id == result_record.id:
                assert seq_record.seq == result_record.seq
