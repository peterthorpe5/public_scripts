#!/usr/bin/env python

"""Test wrapper to /bin/compare_results.py

"""

import os
import subprocess
from pycits.tools import NotExecutableError
from nose.tools import nottest, assert_equal
from pycits.tools import reformat_swarm_cls

INFILE = os.path.join("tests", "test_data", "compare_clusters",
                      "compare_list.txt")
OUTDIR = os.path.join("tests", "test_out_compare")

TARGET = os.path.join("tests", "test_targets", "compare_clusters",
                      "DNAMIX_S95_L001_RESULTS_cd_0.99_sw_1_BC_0.9_V_0.99.txt")


# folder checking
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)


def get_sorted_list(in_file):
    """funct to return a sorted list.
    takes in a file, returns sorted list"""
    out_list = []
    with open(in_file) as fh:
        data = fh.read().split("\n")
        for line in data:
            if not line.strip():
                continue  # if the last line is blank
            if line.startswith("#"):  # dont want comment lines
                continue
        out_list.append(line)
    return sorted(out_list)


def test_compare_results_exec():
    """Run compare_results.py on test data and compare output to
    precomputed target. The default option are the actual
    test data. So we only need to call the program
    """
    outfile = os.path.join(OUTDIR, "compare_tests.txt")
    prog = os.path.join("bin", "compare_results.py")
    temp_s = ["python3",
              prog,
              " -o",
              outfile,
              " --in_list",
              INFILE]
    cmd_s = ' '.join(temp_s)
    pipe = subprocess.run(cmd_s, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
    if not os.path.isfile(outfile):
        sys_exit("outfile not generated: %s" % outfile)
    tests_data = get_sorted_list(TARGET)
    result = get_sorted_list(outfile)
    assert_equal(tests_data, result)
