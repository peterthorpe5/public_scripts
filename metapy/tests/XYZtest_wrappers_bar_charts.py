#!/usr/bin/env python

"""Test wrapper to /bin/draw_bar_chart_of_clusters.py.
The defaults are the actual test data. So it only needs to be called
as METAPY.py
"""

import os
import subprocess
from bin.draw_bar_chart_of_clusters import get_names_from_Seq_db

INFILE = os.path.join("tests", "test_data", "bar_charts",
                      "swarm.outRENAMED_abundance")
INFILE_names = os.path.join("tests", "test_data", "bar_charts",
                            "Names.txt")
INFILE_names_no_abun = os.path.join("tests", "test_data", "bar_charts",
                                    "Names_no_abundance.txt")
OTU_DATABASE_SWARM = os.path.join("tests", "test_data", "bar_charts",
                                  "ITS_db_NOT_confirmed_for_swarm.fasta")

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
    return out_list


def test_get_names():
    """test the get_names_from_Seq_db functio"""
    names, names_abudance_removed = get_names_from_Seq_db(OTU_DATABASE_SWARM)
    test_names = get_sorted_list(INFILE_names)
    test_names_no_abun = get_sorted_list(INFILE_names_no_abun)
    assert_equal(test_names, names)
    assert_equal(test_names_no_abun, names_abudance_removed)


def test_METAPY_exec():
    """Run METAPY.py on test data and compare output to
    precomputed target. The default option are the actual
    test data. So we only need to call the program
    """
    prog = os.path.join("bin", "draw_bar_chart_of_clusters.py")
    temp_s = ["python",
              prog,
              "-i",
              INFILE,
              " --db",
              OTU_DATABASE_SWARM]
    cmd_s = ' '.join(temp_s)
    pipe = subprocess.run(cmd_s, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
    if not os.path.isfile(INFILE + "_barchart.png"):
        sys_exit("outfile not generated: %s" % INFILE + "_barchart.png")
