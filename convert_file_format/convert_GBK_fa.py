#!/usr/bin/env python
# author: Peter Thorpe and Leighton Pritchard. The James Hutton Insitute,Dundee,UK.
# 2017 Feb

# Title:
# script to convert .gbk file to fa.

from Bio import SeqIO

def convert_file(in_file, out_file):
    """function takes in a gbk, writes out a fasta file."""
    SeqIO.convert(in_file, "genbank", out_file, "fasta")
