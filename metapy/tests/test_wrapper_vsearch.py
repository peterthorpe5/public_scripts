#!/usr/bin/env python

"""Tests of wrapper code in pycits. vsearch"""

import os
import shutil
import subprocess

from pycits import vsearch
from pycits.tools import NotExecutableError, reformat_blast6_clusters

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from nose.tools import nottest, assert_equal

# INPUT DATA
INDIR = os.path.join("tests", "test_data", "vsearch")
OUTDIR = os.path.join("tests", "test_out", "vsearch")
TARGETDIR = os.path.join("tests", "test_targets", "vsearch")

TARGET_DEREP = os.path.join(TARGETDIR, "dedup_test_vsearch.fasta")
TARGET_CLUSTER_UC = os.path.join(TARGETDIR, "target.uc")
TARGET_CLUSTER_B6 = os.path.join(TARGETDIR, "target.blast6")
TARGET_R_FORMAT = os.path.join(TARGETDIR, "target.formatR")
TARGET_CLUSTER_FAST_B6 = os.path.join(TARGETDIR, "clusterfast.blast6")
TARGET_CLUSTER_FAST_UC = os.path.join(TARGETDIR,
                                      "target_runfast.clusters.uc")
TARGET_CLUSTER_FAST_MSA = os.path.join(TARGETDIR,
                                      "target.alignedclusters.fasta")
TARGET_CLUSTER_FAST_CENTROIDS = os.path.join(TARGETDIR,
                                             "target.centroids.fasta")
TARGET_CLUSTER_FAST_CONSENSUS = os.path.join(TARGETDIR,
                                             "target.consensus_cls_seq.fasta")

INFILE_DEREP = os.path.join(INDIR, "test_db_for_abundance_checking.fasta")
INFILE_CLUSTER = TARGET_DEREP
DB = os.path.join(INDIR, "vsearch_tests_ITS_db.fasta")
OUTFILE_DEREP = os.path.join(OUTDIR, "vsearch_derep.fasta")
OUTFILE_CLUSTER_UC = os.path.join(OUTDIR, "vsearch_cluster.uc")
OUTFILE_CLUSTER_B6 = os.path.join(OUTDIR, "vsearch_cluster.blast6")
OUTFILE_CLUSTER_FAST_UC = os.path.join(OUTDIR, "vsearch_cluster_fast.uc")
OUTFILE_CLUSTER_FAST_B6 = os.path.join(OUTDIR, "vsearch_cluster_fast.blast6")
OUTFILE_CLUSTER_FAST_MSA = os.path.join(OUTDIR,
                                        "vsearch_cluster_fast_msa.fasta")
OUTFILE_CLUSTER_FAST_CENTROIDS = os.path.join(OUTDIR,
                                              "vsearch_cluster_fast_centroids.fasta")
OUTFILE_CLUSTER_FAST_CONSENSUS = os.path.join(OUTDIR,
                                              "vsearch_cluster_fast_consensus.fasta")

THRESHOLD = 0.96
THREADS = 2

# PARAMETERS
CLUSTER_PARAMS = {'--blast6out': OUTFILE_CLUSTER_B6,
                  '--id': THRESHOLD,
                  '--db': DB,
                  '--threads': THREADS}
CLUSTER_FAST_PARAMS = {"--id": THRESHOLD,
                       "--centroids": OUTFILE_CLUSTER_FAST_CENTROIDS,
                       "--msaout": OUTFILE_CLUSTER_FAST_MSA,
                       "--consout": OUTFILE_CLUSTER_FAST_CONSENSUS,
                       "--db": DB,
                       "--threads": THREADS,
                       "--blast6out": OUTFILE_CLUSTER_FAST_B6}


def setup():
    """Set up test fixtures"""
    try:
        shutil.rmtree(OUTDIR)
    except FileNotFoundError:
        pass
    os.makedirs(OUTDIR, exist_ok=True)


def compare_sorted_lines(file1, file2, skip_comments=True):
    """Compares the sorted line-by-line content of two files.

    If skip_comments is True, then lines starting with '#' are ignored. Lines
    are stripped at both ends for whitespace, and blank lines are discarded.
    """
    with open(file1, 'r') as fh1:
        with open(file2, 'r') as fh2:
            f1_sorted = sorted([l.strip() for l in fh1.readlines() if
                                len(l.strip())])
            f2_sorted = sorted([l.strip() for l in fh2.readlines() if
                                len(l.strip())])
    if skip_comments:
        f1_sorted = [l for l in f1_sorted if not l.startswith('#')]
        f2_sorted = [l for l in f2_sorted if not l.startswith('#')]
    assert_equal(f1_sorted, f2_sorted)


def compare_sorted_fasta(file1, file2):
    """Compares the sorted content of two FASTA files"""
    with open(file1, 'r') as fh1:
        with open(file2, 'r') as fh2:
            fasta1 = list(SeqIO.parse(fh1, 'fasta'))
            fasta2 = list(SeqIO.parse(fh2, 'fasta'))
    fasta1 = sorted([(s.id, s.seq) for s in fasta1])
    fasta2 = sorted([(s.id, s.seq) for s in fasta2])
    assert_equal(fasta1, fasta2)
    
    
def test_vsearch_path():
    """vsearch is in $PATH"""
    vsearch_exe = vsearch.Vsearch("vsearch")


def test_vsearch_exec_notexist():
    """error thrown if vsearch executable does not exist"""
    try:
        vsearch_exe = vsearch.Vsearch(os.path.join(".", "vsearch"))
    except NotExecutableError:
        return True
    else:
        return False


def test_vsearch_notexec():
    """Error thrown if vsearch not executable"""
    try:
        vsearch_exe = vsearch.Vsearch("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_vsearch_derep_cmd():
    """VSEARCH wrapper returns correct dereplication cmd-line"""
    mode = '--derep_fulllength'
    target = ' '.join(['vsearch', mode, INFILE_DEREP,
                       '--output', OUTFILE_DEREP,
                       '--sizeout'])

    vsearch_exe = vsearch.Vsearch("vsearch")
    result = vsearch_exe.run(mode, INFILE_DEREP, OUTFILE_DEREP, dry_run=True)
    assert_equal(result.command, target)


def test_vsearch_derep_exec():
    """VSEARCH dereplicates test data"""
    mode = '--derep_fulllength'
    vsearch_exe = vsearch.Vsearch("vsearch")
    result = vsearch_exe.run(mode, INFILE_DEREP, OUTFILE_DEREP)
    with open(TARGET_DEREP, 'r') as target_fh:
        with open(result.outfilename, 'r') as test_fh:
            assert_equal(target_fh.read(), test_fh.read())

            
def test_vsearch_cluster_cmd():
    """VSEARCH wrapper returns correct clustering cmd-line"""
    mode = '--usearch_global'
    target = ' '.join(['vsearch', mode, INFILE_CLUSTER,
                       '--uc', OUTFILE_CLUSTER_UC,
                       ' '.join(["{0} {1}".format(k, v) for k, v in
                                 sorted(CLUSTER_PARAMS.items())])
    ])

    vsearch_exe = vsearch.Vsearch("vsearch")
    result = vsearch_exe.run(mode, INFILE_CLUSTER, OUTFILE_CLUSTER_UC,
                             CLUSTER_PARAMS, dry_run=True)
    assert_equal(result.command, target)


def test_vsearch_cluster_exec():
    """VSEARCH clusters test data"""
    mode = '--usearch_global'
    vsearch_exe = vsearch.Vsearch("vsearch")
    result = vsearch_exe.run(mode, INFILE_CLUSTER, OUTFILE_CLUSTER_UC,
                             CLUSTER_PARAMS)
    compare_sorted_lines(TARGET_CLUSTER_UC, result.outfile_uc)
    compare_sorted_lines(TARGET_CLUSTER_B6, result.outfile_b6)


def test_vsearch_cluster_fast_cmd():
    """VSEARCH wrapper returns correct cluster_fast cmd-line"""
    mode = '--cluster_fast'
    target = ' '.join(['vsearch', mode, INFILE_CLUSTER,
                       '--uc', OUTFILE_CLUSTER_FAST_UC,
                       ' '.join(["{0} {1}".format(k, v) for k, v in
                                 sorted(CLUSTER_FAST_PARAMS.items())])
    ])

    vsearch_exe = vsearch.Vsearch("vsearch")
    result = vsearch_exe.run(mode, INFILE_CLUSTER,
                             OUTFILE_CLUSTER_FAST_UC,
                             CLUSTER_FAST_PARAMS, dry_run=True)
    assert_equal(result.command, target)


def test_vsearch_cluster_exec():
    """VSEARCH cluster_fast clusters test data"""
    mode = '--cluster_fast'
    vsearch_exe = vsearch.Vsearch("vsearch")
    result = vsearch_exe.run(mode, INFILE_CLUSTER, OUTFILE_CLUSTER_FAST_UC,
                             CLUSTER_FAST_PARAMS)
    compare_sorted_lines(TARGET_CLUSTER_FAST_UC, result.outfile_uc)
    compare_sorted_lines(TARGET_CLUSTER_FAST_B6, result.outfile_b6)
    compare_sorted_fasta(TARGET_CLUSTER_FAST_MSA, result.outfile_msa)
    


###################################################################
# Now to test conversion to another format
@nottest    
def test_convert_vsearch_format():
    """ testing function in tools to convert blast6 format to format
    for R"""
    cat_out = os.path.join(OUTDIR, "db_and_reads.fasta")
    cat_cmd = ' '.join(["cat",
                        DB,
                        INFILE,
                        ">",
                        cat_out])
    pipe = subprocess.run(cat_cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)

    # reformat_blast6_clusters(blast6, db_and_reads, outfile)
    # call the function
    reformat_blast6_clusters(TARGET_BLAST6, cat_out,
                             os.path.join(OUTDIR, "tests.Rformat"))
    # convert to sorted lists
    result_R = get_sorted_list(os.path.join(OUTDIR, "tests.Rformat"))
    target_R = get_sorted_list(TARGET_R_FORMAT)
    assert_equal(result_R, target_R)
