#!/usr/bin/env python

"""Tests pycits wrapper for pear."""

import os
import shutil

from pycits import pear
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

# INPUT DATA
INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_pear")
READS1 = os.path.join(INDIR, "DNAMIX_S95_L001_R1_001.fastq.gz")
READS2 = os.path.join(INDIR, "DNAMIX_S95_L001_R2_001.fastq.gz")
PREFIX = os.path.split(READS1)[-1].split("_R")[0]

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "pear",
                      "pear_test.assembled.fastq")


def test_pear():
    """pear instantiates with cmd-line if pear is in $PATH"""
    pear.Pear("pear")


def test_pear_cmd():
    """pear instantiates, runs and returns correct form of cmd-line"""
    obj = pear.Pear("pear")
    target = ' '.join(["pear", "-f", READS1, "-r", READS2,
                       "--threads", "4", "-o", os.path.join(OUTDIR, PREFIX)])
    assert_equal(obj.run(READS1, READS2, 4, OUTDIR, PREFIX, dry_run=True),
                 target)


def test_pear_exec_notexist():
    """Error thrown if pear executable does not exist"""
    try:
        obj = pear.Pear(os.path.join(".", "pear"))
    except NotExecutableError:
        return True
    else:
        return False


def test_pear_notexec():
    """Error thrown if pear exe not executable"""
    try:
        obj = pear.Pear("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_pear_exec():
    """Run pear on test data and compare output to precomputed target"""
    obj = pear.Pear("pear")
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)
    result = obj.run(READS1, READS2, 4, OUTDIR, PREFIX)
    with open(TARGET, "r") as target_fh:
        with open(result.outfileassembled, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
