#!/usr/bin/env python

"""Tests pycits wrapper for samtools_index.
sort, and IDXSTATS

NOTE TO DO:

consider changing this to PYSAM
https://github.com/pysam-developers/pysam
http://pysam.readthedocs.io/en/latest/"""

import os
import shutil

from pycits import samtools_index
from pycits import samtools_idxstats
from pycits import samtools_sort
from pycits.tools import NotExecutableError
from nose.tools import nottest, assert_equal

THREADS = "1"

# INPUT DATA
INDIR = os.path.join("tests", "test_data", "samtools")
BAM_FILE = os.path.join(INDIR, "Aligned.sortedByCoord.out.bam")

# OUTPUT DATA
OUTDIR = os.path.join("tests", "test_out_smatools")
OUTFILE = os.path.join(OUTDIR, "idxstats")
SORTED_OUTFILE = os.path.join(OUTDIR, "sorted")

# TARGET OUTPUT DATA
TARGET = os.path.join("tests", "test_targets", "samtools",
                      "MisMat_0_star_mappings")


def test_samtools_index():
    """Samtools_index instantiates with cmd-line if Samtools_index
    is in $PATH"""
    samtools_index.Samtools_Index("samtools")


def test_samtools_index_cmd():
    """samtools_index instantiates, runs and returns
    correct form of cmd-line"""
    obj = samtools_index.Samtools_Index("samtools")
    target = ' '.join(["samtools", "index", BAM_FILE])
    stats = obj.run(BAM_FILE)
    assert_equal(stats.command, target)


def test_samtools_index_exec_notexist():
    """Error thrown if samtools_index executable does not exist"""
    try:
        obj = samtools_index.Samtools_Index(os.path.join
                                            (".", "samtools"))
    except NotExecutableError:
        return True
    else:
        return False


def test_samtools_index_notexec():
    """Error thrown if samtools_index exe not executable"""
    try:
        obj = samtools_index.Samtools_Index("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_samtools_index_exec():
    """Run samtools_index on test data and compare output
    to precomputed target"""
    obj = samtools_index.Samtools_Index("samtools")
    result = obj.run(BAM_FILE)
    # test to see if it has produced the bai file.
    if not os.path.isfile(result.bai):
        return False


def test_samtools_idxstats():
    """Samtools_Idxstats instantiates with cmd-line if Samtools_Idxstats
    is in $PATH"""
    samtools_idxstats.Samtools_Idxstats("samtools")


def test_samtools_idxstats_cmd():
    """samtools_idxstats instantiates, runs and returns
    correct form of cmd-line"""
    # checking and making the output folder is staying in for now
    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)
    obj = samtools_idxstats.Samtools_Idxstats("samtools")
    target = ' '.join(["samtools",
                       "idxstats",
                       BAM_FILE,
                       "|",
                       "grep",
                       '-v',
                       "'*'",
                       "|",
                       "sort",
                       "--reverse",
                       "-n",
                       "-k3",
                       ">",
                       OUTFILE])
    stats = obj.run(BAM_FILE, OUTFILE)
    assert_equal(stats.command,
                 target)


def test_samtools_idxstats_exec_notexist():
    """Error thrown if samtools_idxstats executable does not exist"""
    try:
        obj = samtools_idxstats.Samtools_Idxstats(os.path.join
                                                  (".", "samtools"))
    except NotExecutableError:
        return True
    else:
        return False


def test_samtools_idxstats_notexec():
    """Error thrown if samtools_idxstats exe not executable"""
    try:
        obj = samtools_idxstats.Samtools_Idxstats("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_samtools_idxstats_exec():
    """Run samtools_idxstats on test data and compare output
    to precomputed target"""
    obj = samtools_idxstats.Samtools_Idxstats("samtools")
    result = obj.run(BAM_FILE, OUTFILE)
    with open(TARGET, "r") as target_fh:
        with open(result.cov, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())


def test_samtools_sort_cmd():
    """samtools_sort instantiates, runs and returns
    correct form of cmd-line"""
    obj = samtools_sort.Samtools_Sort("samtools")
    target = ' '.join(["samtools",
                       "sort",
                       "-@",
                       THREADS,
                       BAM_FILE,
                       '-o', SORTED_OUTFILE])
    results = obj.run(BAM_FILE, SORTED_OUTFILE, THREADS)
    assert_equal(results.command,
                 target)


def test_samtools_index_exec():
    """Run samtools_index on test data and compare output
    to precomputed target"""
    obj = samtools_sort.Samtools_Sort("samtools")
    result = obj.run(BAM_FILE, SORTED_OUTFILE, THREADS)
    # test to see if it has produced the bai file.
    if not os.path.isfile(result.sorted_bam):
        return False
