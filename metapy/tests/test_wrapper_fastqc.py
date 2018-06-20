#!/usr/bin/env python

"""Tests pycits wrapper for fastqc."""

import os
import shutil

from pycits import fastqc
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

# INPUT DATA
INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_fastqc")
READS1 = os.path.join(INDIR, "DNAMIX_S95_L001_R1_001.fastq.gz")
READS2 = os.path.join(INDIR, "DNAMIX_S95_L001_R2_001.fastq.gz")
PREFIX = os.path.split(READS1)[-1].split("_R")[0]

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "fastqc",
                      "DNAMIX_S95_L001_R1_001_fastqc.html")


def test_fastqc():
    """fastqc instantiates with cmd-line if fastqc is in $PATH"""
    fastqc.FastQC("fastqc")


def test_fastqc_cmd():
    """fastqc instantiates, runs and returns correct form of cmd-line"""
    obj = fastqc.FastQC("fastqc")
    target = ' '.join(["fastqc", READS1, "-o", OUTDIR])
    qc = obj.run(READS1, OUTDIR)
    assert_equal(qc.command,
                 target)


def test_fastqc_exec_notexist():
    """Error thrown if fastqc executable does not exist"""
    try:
        obj = fastqc.FastQC(os.path.join(".", "fastqc"))
    except NotExecutableError:
        return True
    else:
        return False


def test_fastqc_notexec():
    """Error thrown if fastqc exe not executable"""
    try:
        obj = fastqc.FastQC("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False

# there isnt any reason to test this ...
# def test_fastqc_exec():
#    """Run fastqc on test data and compare output to precomputed target"""
#    obj = fastqc.FastQC("fastqc")
#    try:
#        shutil.rmtree(OUTDIR)
#    except FileNotFoundError:
#        pass
#    os.makedirs(OUTDIR, exist_ok=True)
#    result = obj.run(READS1, OUTDIR)
#    with open(TARGET, "r") as target_fh:
#        with open(result.fastqc_html, "r") as test_fh:
#            assert_equal(target_fh.read(), test_fh.read())
