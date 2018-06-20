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
INDIR = os.path.join("tests", "test_data", "bowtie2")
OUTDIR = os.path.join("tests", "test_out", "bowtie2")
TARGET = os.path.join("tests", "test_targets", "bowtie2",
                      "bowtie2_target_mapping.sam")
READS = os.path.join(INDIR, "bowtie_test.assembled.fasta")
FA_INDEX = os.path.join("tests", "test_data", "bowtie2", "fasta_index")
THREADS = "1"


# Create output directory tree
def setup():
    """Set up test fixtures"""
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)


def test_bowtie2path():
    """bowtie2-map executable is in $PATH"""
    bowtie_map.Bowtie2_Map("bowtie2")


def test_bowtie2_cmd():
    """bowtie2-map returns correct form of cmd-line"""
    bt2_map = bowtie_map.Bowtie2_Map("bowtie2")
    outfilename = os.path.join(OUTDIR,
                               "_vs_".join([os.path.split(os.path.splitext(READS)[0])[-1],
                                            os.path.split(os.path.splitext(FA_INDEX)[0])[-1]]) + ".sam")
    target = ' '.join(["bowtie2",
                       "--very-sensitive",
                       "--no-unal",
                       "-p", THREADS,
                       "-x", FA_INDEX,
                       "-f", READS,
                       "-S", outfilename])
    result = bt2_map.run(READS.split(".gz")[0], FA_INDEX,
                         outfilename, THREADS, fasta=True, dry_run=True)
    assert_equal(result.command, target)


def test_bowtie2_exec_notexist():
    """Error thrown if bowtie2 executable does not exist"""
    try:
        obj = bowtie_map.Bowtie2_Map(os.path.join(".", "bowtie2"))
    except NotExecutableError:
        return True
    else:
        return False


def test_bowtie2_notexec():
    """Error thrown if bowtie2 not executable"""
    try:
        obj = bowtie_map.Bowtie2_Map("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_bowtie2_build_exec():
    """bowtie2 maps reads correctly with FASTA test data"""
    bt2_map = bowtie_map.Bowtie2_Map("bowtie2")
    outfilename = os.path.join(OUTDIR,
                               "_vs_".join([os.path.split(os.path.splitext(READS)[0])[-1],
                                            os.path.split(os.path.splitext(FA_INDEX)[0])[-1]]) +
                               ".sam")
    result = bt2_map.run(READS, FA_INDEX, outfilename, THREADS, fasta=True)
    # Test for equality of output and target MD5 files
    # To establish equality, we need to ignore the @PG line that describes the
    # path to bowtie2, as this may differ between systems
    with open(os.path.join(outfilename), "r") as outfh:
        with open(TARGET, "r") as tgtfh:
            tgtdata = '\n'.join([l for l in tgtfh.readlines() if
                                 not l.startswith('@PG')])
            outdata = '\n'.join([l for l in outfh.readlines() if
                                 not l.startswith('@PG')])
            assert_equal(tgtdata, outdata)
