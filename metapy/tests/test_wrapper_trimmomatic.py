#!/usr/bin/env python

"""Tests pycits wrapper for Trimmomatic"""

import hashlib
import os
import shutil
import subprocess
import sys

from pycits import trimmomatic
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal

# INPUT DATA LOCATIONS
INDIR = os.path.join("tests", "test_data")
OUTDIR = os.path.join("tests", "test_out_trimmomatic")
READS1 = os.path.join(INDIR, "DNAMIX_S95_L001_R1_001.fastq.gz")
READS2 = os.path.join(INDIR, "DNAMIX_S95_L001_R2_001.fastq.gz")
PREFIX = os.path.split(READS1)[-1].split("_R")[0]
PHREDSCORE = "phred33"
ADAPTERS = os.path.join(INDIR, "adapters", "TruSeq3-PE.fa")
OUTFILES = [os.path.join(OUTDIR, PREFIX + suffix) for suffix in
            ("_paired_R1.fq.gz", "_unpaired_R1.fq.gz",
             "_paired_R2.fq.gz", "_unpaired_R2.fq.gz")]

# TARGET OUTPUT
TARGETDIR = os.path.join("tests", "test_targets", "trimmomatic")
TARGETS = [os.path.join(TARGETDIR, PREFIX + suffix) for suffix in
           ("_paired_R1.fq.gz", "_unpaired_R1.fq.gz",
            "_paired_R2.fq.gz", "_unpaired_R2.fq.gz")]


def test_trimmomatic():
    """Trimmomatic instantiates with cmd-line and trimmomatic is in $PATH

    NOTE: there are several ways to install Trimmomatic. We assume that
    trimmomatic is a valid command in the $PATH. This is how homebrew and
    apt-get install trimmomatic. If a 'traditional' .jar installation is
    all that is available, then the command line is
    "java -jar <PATH TO TRIMMOMATIC>" and should be substituted
    """
    trimmomatic.Trimmomatic("trimmomatic")


def test_trimmomatic_cmd():
    """Trimmomatic instantiates and returns correct form of cmd-line"""
    trim = trimmomatic.Trimmomatic("trimmomatic")

    parameters = trimmomatic.Parameters(threads=4)
    steps = trimmomatic.Steps(ILLUMINACLIP="{0}:2:30:10".format(ADAPTERS))
    target = ' '.join(["trimmomatic", "PE", "-phred33", "-threads 4",
                       READS1, READS2, *OUTFILES,
                       "ILLUMINACLIP:{0}:2:30:10".format(ADAPTERS),
                       "LEADING:3", "HEADCROP:0", "TRAILING:3",
                       "MINLEN:70", "SLIDINGWINDOW:4:25"])

    assert_equal(trim.run(READS1, READS2, OUTDIR, PREFIX, PHREDSCORE,
                          parameters, steps, dry_run=True),
                 target)


def test_trimmomatic_exec():
    """Run trimmomatic on test data"""
    trim = trimmomatic.Trimmomatic("trimmomatic")
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)

    parameters = trimmomatic.Parameters(threads=4)
    steps = trimmomatic.Steps(ILLUMINACLIP="{0}:2:30:10".format(ADAPTERS))
    results = trim.run(READS1, READS2, OUTDIR, PREFIX, PHREDSCORE,
                       parameters, steps)

    for outfname, targetfname in zip((results.outfileR1paired,
                                      results.outfileR1unpaired,
                                      results.outfileR2paired,
                                      results.outfileR2unpaired), TARGETS):
        with open(outfname, 'rb') as ofh:
            with open(targetfname, 'rb') as tfh:
                ohash = hashlib.md5(ofh.read()).digest()
                thash = hashlib.md5(tfh.read()).digest()
                assert_equal(ohash, thash)
