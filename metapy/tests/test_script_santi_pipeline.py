#!/usr/bin/env python3

"""Tests of run_santi_pipeline.py script."""

import os
import shutil

from argparse import Namespace
from nose.tools import nottest, assert_equal, with_setup

from pycits.scripts import santi_pipeline

NAMESPACE = Namespace(blastclust='blastclust',
                      convert_format='convert_format',
                      fastqc='fastqc', indirname='data/reads',
                      join_paired_ends='join_paired_ends.py',
                      logfile='nosetest_santi_script_test.log',
                      muscle='muscle',
                      outdirname='nosetest_santi_script_test',
                      pick_closed_reference_otus='pick_closed_reference_otus.py',
                      pick_otus='pick_otus.py',
                      prefix='DNAMIX_S95_L001',
                      reference_fasta='data/database.fasta',
                      threads=4,
                      trim_quality='trim_quality',
                      verbose=True)
DATADIR = os.path.join("tests", "test_data", "santi_pipeline")
OUTDIR = os.path.join("tests", "test_out", "santi_pipeline")
TARGETDIR = os.path.join("tests", "test_targets", "santi_pipeline")


# The setup_outdir() function decorates functions by creating the appropriate
# output directory tree
def setup_outdir():
    """Set up test fixtures"""
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)


def test_complete_script_notravis():
    """complete script runs"""
    if os.path.isdir(NAMESPACE.outdirname):
        shutil.rmtree(NAMESPACE.outdirname)
    santi_pipeline.run_pycits_main(NAMESPACE)


@with_setup(setup_outdir)
def test_biom_to_tsv():
    """biom format to TSV conversion"""
    infname = os.path.join(DATADIR, "otu_table.biom")
    outdir = os.path.join(OUTDIR)
    targetfname = os.path.join(TARGETDIR, "otu_table.tsv")
    outfname = santi_pipeline.biom_to_tsv(infname, outdir)
    with open(targetfname, "r") as target_fh:
        with open(outfname, "r") as output_fh:
            assert_equal(target_fh.read(), output_fh.read())


def test_logger_creation():
    """santi_pipeline log creation"""
    logger = santi_pipeline.construct_logger(NAMESPACE, header=False)


def test_outdir_creation():
    """santi_pipeline creates outdir"""
    outdirname = os.path.join(OUTDIR, "test_outdir")
    santi_pipeline.create_output_directory(outdirname)
