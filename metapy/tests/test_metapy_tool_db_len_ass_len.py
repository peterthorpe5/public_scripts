#!/usr/bin/env python3

"""Test program for a function in the module metapy_tools:
db_len_assembled_len_ok
"""

import os
import shutil

from nose.tools import nottest, assert_equal
from pycits.metapy_tools import db_len_assembled_len_ok
from pycits.tools import NotExecutableError

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data", "metapy_tools")
DB = os.path.join("data", "ITS_database_NOT_confirmed_correct" +
                  "_last14bases_removed.fasta")
OUTDIR = os.path.join("tests", "test_out", "metapy_tools")
INFILE_TOO_SHORT = os.path.join(INDIR,
                                "assembled_fa_too_short.fasta")
INFILE_TOO_LONG = os.path.join(INDIR,
                               "assembled_fa_too_long.fasta")
INFILE_GOOD = os.path.join(INDIR,
                           "assembled_fa_good.fasta")

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


def test_assem_db_seq_len_short():
    """Run db_len_assembled_len_ok too SHORT from metapy_tools
    """
    error, msg = db_len_assembled_len_ok(DB,
                                         INFILE_TOO_SHORT)
    errtype, db_mean, db_sd, assemb_mean, assemb_sd = msg.split("\t")
    assert error == "fail"
    # msg : -       196.139 18.997  21.636  8.391
    assert errtype == "-"
    assert str(db_mean) == "196.139"
    assert str(db_sd) == "18.997"
    assert str(assemb_mean) == "21.636"
    assert str(assemb_sd) == "8.391"


def test_assem_db_seq_len_long():
    """Run db_len_assem_len_ok too LONG from metapy_tools
    """
    error, msg = db_len_assembled_len_ok(DB,
                                         INFILE_TOO_LONG)
    errtype, db_mean, db_sd, assemb_mean, assemb_sd = msg.split("\t")
    assert error == "fail"
    # msg : "+       196.139 18.997  1150.800        141.426"
    assert errtype == "+"
    assert str(db_mean) == "196.139"
    assert str(db_sd) == "18.997"
    assert str(assemb_mean) == "1150.800"
    assert str(assemb_sd) == "141.426"


def test_assem_db_seq_len_ok():
    """Run db_len_assem_len_ok OK from metapy_tools
    """
    error, msg = db_len_assembled_len_ok(DB,
                                         INFILE_GOOD,
                                         10)
    # ok   196.139 18.997  182.684 19.263
    assert error == "ok"
    errtype, db_mean, db_sd, assemb_mean, assemb_sd = msg.split("\t")
    assert errtype == "ok"
    assert str(db_mean) == "196.139"
    assert str(db_sd) == "18.997"
    assert str(assemb_mean) == "182.684"
    assert str(assemb_sd) == "19.263"


def test_assem_db_seq_len_ok_sd_0():
    """Run db_len_assem_len 0 sd threshold
    """
    error, msg = db_len_assembled_len_ok(DB,
                                         INFILE_GOOD,
                                         0)
    assert error == "fail"
