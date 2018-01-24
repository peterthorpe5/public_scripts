#!/usr/bin/env python
# title write a shelll to run metapy on mutli fastaq files.

import os
import shutil
import subprocess
import sys
from sys import stdin,argv

wanted_name = "'" + argv[1] + "*'"
file_names = open("all_results_files.txt", "r")
command = " ".join(["echo",
                        "'#filename species\tNumber_of_reads'"
                         ">",
                         argv[1] + "identified.out"])
pipe = subprocess.run(command, shell=True,
                      stdout=subprocess.PIPE,
                      stderr=subprocess.PIPE)
for filename in file_names:
    command = " ".join(["grep",
                        "-H", 
                        "--with-filename",
                         wanted_name,
                         filename.rstrip(),
                        "|",
                        "grep",
                        "-v",
                        "'cd_hit'",
                        "|",
                        "grep",
                        "-v",
                        "'blastclust'",
                        "|",
                        "grep",
                        "-v",
                        "'vserach_clu_fasta'",
                         ">>",
                         argv[1] + "identified.out"])
    # print(command)
    pipe = subprocess.run(command, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    # print(pipe.stdout, pipe.stderr)
