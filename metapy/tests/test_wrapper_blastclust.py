#!/usr/bin/env python

"""Tests pycits wrapper for blastclust.
Also tests the conversion to another format"""

import os
import shutil

from pycits import blast
from pycits.tools import NotExecutableError
from pycits.tools import reformat_swarm_cls

from nose.tools import nottest, assert_equal

THREADS = 2

# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data", "blastclust")
OUTDIR = os.path.join("tests", "test_out", "blastclust")

# TARGET OUTPUT
TARGET = os.path.join("tests", "test_targets", "blastclust",
                      "target_trimmed.fasta.blastclust99.lst")

TARGET_FOR_R = os.path.join("tests", "test_targets", "blastclust",
                            "trimmed.fasta.blastclust99.Rout")


# Create output directory tree
def setup():
    """Set up test fixtures"""
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)


def test_blastclust():
    """blastclust is in $PATH"""
    blast.Blastclust("blastclust")


def test_blastclust_cmd():
    """Blastclust instantiates and returns correct form of cmd-line"""
    bc = blast.Blastclust("blastclust")
    target = ' '.join(["blastclust", "-L 0.90", "-S 99", "-a 4", "-p F",
                       "-i trimmed.fasta",
                       "-o", os.path.join("tests", "test_out", "blastclust",
                                          "trimmed.fasta.blastclust99.lst")])
    result = bc.run("trimmed.fasta", OUTDIR, 4, dry_run=True)
    assert_equal(result.command, target)


def test_blastclust_exec_notexist():
    """Error thrown if executable does not exist"""
    try:
        bc = blast.Blastclust(os.path.join(".", "blastclust"))
    except NotExecutableError:
        return True
    else:
        return False


def test_blastclust_notexec():
    """Error thrown if blastclust exe not executable"""
    try:
        bc = blast.Blastclust("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_blastclust_exec():
    """Run blastclust on test data and compare output to precomputed target"""
    bc = blast.Blastclust("blastclust")
    result = bc.run(os.path.join(INDIR, "trimmed.fasta"), OUTDIR, THREADS)
    with open(TARGET, "r") as target_fh:
        with open(result.outfilename, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
