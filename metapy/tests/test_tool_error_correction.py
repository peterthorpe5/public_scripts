#!/usr/bin/env python

"""Tests of wrapper code in pycits. spades.py for
error correction
"""

import os
import shutil
import gzip
from pycits import error_correction
from pycits.tools import NotExecutableError
from nose.tools import nottest, assert_equal
from Bio.SeqIO.QualityIO import FastqGeneralIterator


# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_EC")
READS1 = os.path.join(INDIR, "DNAMIX_S95_L001_R1_001.fastq.gz")
READS2 = os.path.join(INDIR, "DNAMIX_S95_L001_R2_001.fastq.gz")
PREFIX = "DNAMIX_S95_L001_"
# TARGET OUTPUT DATA
TARGET_L = os.path.join("tests", "test_targets", "error_correction",
                        "DNAMIX_S95_L001_R1_001.fastq.00.0_0.cor.fastq.gz")
TARGET_R = os.path.join("tests", "test_targets", "error_correction",
                        "DNAMIX_S95_L001_R2_001.fastq.00.0_0.cor.fastq.gz")


def sort_fq_output(fq_file):
    """function to sort the output of the EC assembled data
    . Normally this is returned in an unordered manner. Needs
    sorting for tests.
    returned as a sorted(set)"""
    # open the file
    try:
        # the files should be .gz, py lib to open
        in_file = gzip.open(fq_file, mode='rt', compresslevel=9,
                            encoding=None, errors=None,
                            newline=None)
    except ValueError:
        fq_file = fq_file.split(".gz")[0]
        in_file = open(fq_file)
    fq_set = set([])
    for (title, seq, qual) in (FastqGeneralIterator(in_file)):
        read_info = "%s\n+\n%s\n" % (seq, qual)
        fq_set.add(title + read_info)
    in_file.close()
    return sorted(fq_set)


def test_Error_Correction():
    """spades.py instantiates with cmd-line if spades.py is in $PATH"""
    ec = error_correction.Error_Correction("spades.py")


def test_Error_Correction_cmd():
    """spades.py instantiates, runs and returns correct form of cmd-line"""
    obj = error_correction.Error_Correction("spades.py")
    target = ' '.join(["spades.py",
                       "-1",
                       READS1,
                       "-2",
                       READS2,
                       "--only-error-correction",
                       "--threads", "4",
                       "-o", OUTDIR])
    assert_equal(obj.run(READS1, READS2, 4,
                         os.path.join(OUTDIR),
                         dry_run=True), target)


def test_spade_py_exec_notexist():
    """Error thrown if spade.py executable does not exist"""
    try:
        ec = error_correction.Error_Correction(os.path.join(".", "spades.py"))
    except NotExecutableError:
        return True
    else:
        return False


def test_EC_notexec():
    """Error thrown if EC exe not executable"""
    try:
        obj = error_correction.Error_Correction("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_error_correction_exec():
    """Run spades.py on test data and compare output to precomputed target
    """
    obj = error_correction.Error_Correction("spades.py")
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)
    result = obj.run(READS1, READS2, 4, OUTDIR, dry_run=False)
    # print ("results:", result.right_read_correct, "\n\n")
    # call function to sort the data. returned as a sorted(set)
    # the output is not sorted, so a direct comparason of the file
    # does not work. So have to sort first.
    # print ("now sorting TARGET_L", TARGET_L, "\n\n")
    targed_data_sorted_left = sort_fq_output(TARGET_L)
    # print ("now sorting TARGET_R", TARGET_R, "\n\n")
    targed_data_sorted_right = sort_fq_output(TARGET_R)
    # Left_read_correct - name from named tuple and originally
    # print ("sorting result.Left_read_correct", result.Left_read_correct)
    test_data_sorted_L = sort_fq_output(result.Left_read_correct)
    test_data_sorted_R = sort_fq_output(result.right_read_correct)
    # test if they are equal
    assert_equal(targed_data_sorted_left, test_data_sorted_L)
    assert_equal(targed_data_sorted_right, test_data_sorted_R)
