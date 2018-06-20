#!/usr/bin/env python

"""Tests pycits wrapper for muscle."""

import os
import shutil

from pycits import muscle
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal, with_setup

# INPUT DATA
INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out", "muscle")
INFILE = os.path.join(INDIR, "dedup_test.fasta")
OUTFILE = os.path.join(OUTDIR, "dedup_test_aln.fasta")

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "muscle",
                      "muscle_tests.fasta")


# The setup() function runs before all the tests (teardown() runs after they
# are complete).
def setup():
    """Set up test fixtures"""
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)


def test_muscle_path():
    """muscle executable is in $PATH"""
    muscle_aln = muscle.Muscle("muscle")


def test_muscle_exec_notexist():
    """error thrown if muscle executable does not exist"""
    try:
        bc = muscle.Muscle(os.path.join(".", "muscle"))
    except NotExecutableError:
        return True
    else:
        return False


def test_muscle_cmd():
    """muscle wrapper returns correct form of cmd-line"""
    muscle_aln = muscle.Muscle("muscle")
    muscle_result = muscle_aln.run(INFILE, OUTFILE, dry_run=True)
    target = ' '.join(['muscle',
                       '-in', INFILE,
                       '-out', OUTFILE])
    assert_equal(muscle_result.command, target)


def test_muscle_exec():
    """muscle aligns test data"""
    muscle_aln = muscle.Muscle("muscle")
    result = muscle_aln.run(INFILE, OUTFILE)
    with open(TARGET, "r") as target_fh:
        with open(result.outfile, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())


def test_muscle_notexec():
    """error thrown if muscle exe not executable"""
    try:
        bc = muscle.Muscle("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False
