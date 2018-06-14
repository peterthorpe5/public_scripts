#!/usr/bin/env python
# coding: utf-8
# Title:
# script to reformat fastq to fasta format
# Why: after pear/ flash have assembled the reads.
# Need to convert from fastq to fa
# author: Peter Thorpe and Leighton Pritchard
# November 2016. The James Hutton Insitute, Dundee, UK.


from Bio import SeqIO

#os imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser

def convert_file(in_file, out_file):
    SeqIO.convert(in_file, "fastq", out_file, "fasta")

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

$ convert_fq_to_fa -i in.fastq -o outfile.fasta

requires bippython
"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="in_file", default=None,
                  help="fastq file")
parser.add_option("-o", "--output", dest="out_file", default=None,
                  help="Output filename, fasta",
                  metavar="FILE")


(options, args) = parser.parse_args()

in_file = options.in_file
out_file = options.out_file

(options, args) = parser.parse_args()
# run the program
convert_file(in_file, out_file)
