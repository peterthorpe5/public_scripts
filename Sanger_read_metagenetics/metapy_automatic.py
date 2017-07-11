#!/usr/bin/env python
#
# metapy_sanger_read.py -  front end script to run over everything in a
# folder
#
# Code for script to identify OTUs from metabarcoding reads.
# runs multiple clustering programs
# THIS IS FOR SANGER READS
#
# (c) The James Hutton Institute 2016-2017
# Author: Leighton Pritchard, Peter Thorpe

import os
import subprocess

cwd = os.getcwd()

for filename in os.listdir(".") :
    if not filename.endswith(".ab1"):
        continue
    python_cmd = " ".join([os.path.join("$HOME",
                                        "public_scripts",
                                        "Sanger_read_metagenetics",
                                        "metapy_sanger_read.py"),
                           "-a",
                           filename])
    print("running command: %s" % python_cmd)
    pipe = subprocess.run(python_cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    print("finished %s" % filename)

print("finsihed")
