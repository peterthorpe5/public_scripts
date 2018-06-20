#!/usr/bin/env python

"""Tests pycits wrapper for bowtie2.

Bowtie does not have native support for gzipped files
"""

import hashlib
import os
import shutil

from pycits import bowtie_build, bowtie_map
from pycits.tools import NotExecutableError
from nose.tools import nottest, assert_equal


import subprocess

# INPUT DATA
OUTDIR = os.path.join("tests", "test_out", "bowtie2-build")
TARGETDIR = os.path.join("tests", "test_targets", "bowtie2-build")
DATABASE = os.path.join("tests", "test_data", "bowtie2-build",
                        "database_bowtie_test.fasta")

FA_INDEX = os.path.join(OUTDIR, "fasta_index")
THREADS = "1"


# Create output directory tree
def setup():
    """Set up test fixtures"""
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)


# test for indexing first
def test_bowtie2_build_path():
    """bowtie2-build executable is in $PATH"""
    bowtie_build.Bowtie2_Build("bowtie2-build")


def test_bowtie2_build_exec_notexist():
    """Error thrown if bowtie_build executable does not exist"""
    try:
        obj = bowtie_build.Bowtie2_Build(os.path.join(".", "bowtie2-build"))
    except NotExecutableError:
        return True
    else:
        return False


def test_bowtie2_build_cmd():
    """bowtie2-build returns correct form of cmd-line"""
    bowtie2_idx = bowtie_build.Bowtie2_Build("bowtie2-build")
    target = ' '.join(["bowtie2-build",
                       "--quiet",
                       "-f",
                       DATABASE,
                       FA_INDEX])
    result = bowtie2_idx.run(DATABASE, FA_INDEX, dry_run=True)
    assert_equal(result.command, target)


def test_bowtie2_build_notexec():
    """Error thrown if bowtie2-build not executable"""
    try:
        obj = bowtie_build.Bowtie2_Build("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_bowtie2_build_exec():
    """bowtie2-build indexes the file correctly"""
    bowtie2_idx = bowtie_build.Bowtie2_Build("bowtie2-build")
    result = bowtie2_idx.run(DATABASE, FA_INDEX)
    # Test for equality of output and target MD5 hashes
    for fname in os.listdir(OUTDIR):
        with open(os.path.join(OUTDIR, fname), "rb") as outfh:
            outhash = hashlib.md5()
            outhash.update(outfh.read())
        with open(os.path.join(TARGETDIR, fname), "rb") as tgtfh:
            tgthash = hashlib.md5()
            tgthash.update(tgtfh.read())
        assert_equal(tgthash.digest(), outhash.digest())
