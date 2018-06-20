#!/usr/bin/env python

"""Tests of wrapper code in pycits. cd_hit_est"""

import os
import shutil

from pycits import cd_hit
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

# INPUT DATA LOCATION
INDIR = os.path.join("tests", "test_data", "cd_hit")
OUTDIR = os.path.join("tests", "test_out_cd_hit")
INFILE = os.path.join(INDIR, "assembled_reads_and_OTU_db.fasta")
OUTFILE = os.path.join(OUTDIR, "cd_hit_tests.out.clstr")
PREFIX = "test_run"

# TARGET OUTPUT
TARGET = os.path.join("tests", "test_targets", "cd_hit",
                      "cd_hit_tests.out.clstr")


def test_cd_hit():
    """cd_hit instantiates with cmd-line if cd-hit is in $PATH"""
    cluster = cd_hit.Cd_hit("cd-hit-est")


def test_cd_hit_cmd():
    """cd_hit instantiates and returns correct form of cmd-line"""
    cluster = cd_hit.Cd_hit("cd-hit-est")
    target = ' '.join(["cd-hit-est",
                       "-i", INFILE,
                       "-o", os.path.join(OUTDIR, PREFIX),
                       "-T", "4",
                       "-M", "0",
                       "-c", "0.99",
                       "-d", "500"])
    assert_equal(cluster.run(INFILE, "4", "0.99", OUTDIR, PREFIX,
                             dry_run=True), target)


def test_cd_hit_exec_notexist():
    """Error thrown if CD-Hit executable does not exist"""
    try:
        cluster = cd_hit.Cd_hit(os.path.join(".", "cd-hit-est"))
    except NotExecutableError:
        return True
    else:
        return False


def test_cd_hit_notexec():
    """Error thrown if cd_hit not executable"""
    try:
        cluster = cd_hit.Cd_hit("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_cd_hit_exec():
    """Run cd_hit on test data

    TODO: finer option could be passed??
    """
    cluster = cd_hit.Cd_hit("cd-hit-est")
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)

    # results are defined in the class:
    # factory class for cd_hit class returned values
    # Results = namedtuple("Results", "command fastaout clusters " +
    # "stdout stderr")
    # fasta_in, threads, threshold, outdir, prefix
    result = cluster.run(INFILE, "4", "0.99", OUTDIR, PREFIX, dry_run=False)
    # use the named tuple to get the cluster results file
    cd_hit_clusters = result.clusters

    with open(TARGET, "rt") as target_fh:
        with open(cd_hit_clusters, "r") as test_fh:
            assert_equal(target_fh.read(), test_fh.read())
