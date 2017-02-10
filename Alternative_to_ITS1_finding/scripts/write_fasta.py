#!/usr/bin/env python
# author: Peter Thorpe and Leighton Pritchard. The James Hutton Insitute,Dundee,UK.
# 2017 Feb

# Title: simply to write the result out as a fasta file

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
   

def write_fasta(filename, outfile):
    "this function re-write a file as a fasta file"
    # outfile is already opened
    for seq_record in SeqIO.parse(filename, "fasta"):
        SeqIO.write(seq_record, outfile, "fasta")
    # do not close the outfile here!!

