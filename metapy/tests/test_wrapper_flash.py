#!/usr/bin/env python

"""Tests of wrapper code in pycits: Flash"""

import os
import shutil

from pycits import flash
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

from Bio.SeqIO.QualityIO import FastqGeneralIterator

# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_flash")
READS1 = os.path.join(INDIR, "DNAMIX_S95_L001_R1_001.fastq.gz")
READS2 = os.path.join(INDIR, "DNAMIX_S95_L001_R2_001.fastq.gz")
PREFIX = os.path.split(READS1)[-1].split("_R")[0]

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "flash",
                      "flash_test.assembled.fastq")

# folder checking
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)


def sort_flash_output(flash_assmebled_file):
    """function to sort the output of the flash assembled data
    . Normally this is returned in an unordered manner. Needs
    sorting for tests.
    returned as a sorted(set)"""
    # open the file
    in_file = open(flash_assmebled_file)
    fq_set = set([])
    for (title, seq, qual) in (FastqGeneralIterator(in_file)):
        read_info = "%s\n+\n%s\n" % (seq, qual)
        fq_set.add(title + read_info)
    in_file.close()
    return sorted(fq_set)


def test_flash():
    """flash instantiates with cmd-line if flash is in $PATH"""
    assemble = flash.Flash("flash")


def test_flash_cmd():
    """flash instantiates, runs and returns correct form of cmd-line"""
    obj = flash.Flash("flash")
    target = ' '.join(["flash", "--max-overlap", "250",
                       "--threads", "4", "-o", os.path.join(OUTDIR, PREFIX),
                       READS1, READS2])
    assert_equal(obj.run(READS1, READS2, 4, OUTDIR, PREFIX, dry_run=True),
                 target)


def test_flash_exec_notexist():
    """Error thrown if flash executable does not exist"""
    try:
        obj = flash.Flash(os.path.join(".", "flash"))
    except NotExecutableError:
        return True
    else:
        return False


def test_flash_notexec():
    """Error thrown if flash exe not executable"""
    try:
        obj = flash.Flash("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_flash_exec():
    """Run flash on test data and compare output to precomputed target
    this often does not produce the same result twice. So we may have
    to disable this test"""
    obj = flash.Flash("flash")
    result = obj.run(READS1, READS2, 4, OUTDIR, PREFIX, dry_run=False)
    # call function to sort the data. returned as a sorted(set)
    # the flash output is not sorted, so a direct comparason of the file
    # does not work. So have to sort first.
    targed_data_sorted = sort_flash_output(TARGET)
    # outfileextended - name from named tuple and originally
    # from build_command
    test_data_sorted = sort_flash_output(result.outfileextended)
    assert_equal(targed_data_sorted, test_data_sorted)
