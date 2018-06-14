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
# Author: Peter Thorpe, Leighton Pritchard

import os
import subprocess
import sys

if sys.version_info[:2] != (3, 5):
    # e.g. sys.version_info(major=3, minor=5, micro=2,
    # releaselevel='final', serial=0)
    # break the program
    if sys.version_info[:2] != (2, 7):
        print ("currently using:", sys.version_info,
               "  version of python")
        raise ImportError("Python 3.5 or 2.7 is required for " +
                          "metapy_sanger_read.py")
        print ("did you activate the virtual environment?")
        sys.exit(1)

VERSION = """Pycits/ metapy classify OTU using Sanger ab1 files: v1.0.0.
This script will run the method of all .abi or .ab1 files in a folder. 
This come with a Phytophora database, which are all the enteries in 
NCBI as of Jan 2017. If you want your own database you will have to run 
\nmetapy_sanger_read.py -d your_datase.fasta\n"""
if "--version" in sys.argv:
    print(VERSION)
    sys.exit(1)

cwd = os.getcwd()

for filename in os.listdir(".") :
    if not filename.endswith(".ab1"):
        continue
    python_cmd = " ".join([os.path.join("$HOME",
                                        "public_scripts",
                                        "Sanger_read_metagenetics",
                                        "metapy_sanger_read.py"),
                           "-a",
                           filename,
                           "--thread",
                           "16"])
                           #"-d", "~/misc_python/THAPBI/THAPBI-pycits/data/Phytophora_ITS_database_v0.004.fasta"])

    #"-d", "~/misc_python/THAPBI/THAPBI-pycits/data/Phytophora_ITS_database_v0.004.fasta"
    #default is nt database.
    print("running command: %s" % python_cmd)
    pipe = subprocess.run(python_cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    print("finished %s" % filename)

for filename in os.listdir(".") :
    if not filename.endswith(".abi"):
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
