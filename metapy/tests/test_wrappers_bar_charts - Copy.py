#!/usr/bin/env python

"""Test wrapper to /bin/draw_bar_chart_of_clusters.py.
The defaults are the actual test data. So it only needs to be called
as METAPY.py
"""

import os
import subprocess

INFILE = os.path.join("tests", "test_data", "bar_charts",
                      "swarm.outRENAMED_abundance")

OTU_DATABASE_SWARM = os.path.join("tests", "test_data", "bar_charts",
                                  "ITS_db_NOT_confirmed_for_swarm.fasta")


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
