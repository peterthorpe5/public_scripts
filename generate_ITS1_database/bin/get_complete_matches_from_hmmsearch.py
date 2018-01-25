#!/usr/bin/env python

# title: script to get the complete matches from hmmsearch
# which do not get put in the dmtblout file
# e.g. >> 9_Phytophthora_hydropathica_HG934150  
#  [No individual domains that satisfy reporting thresholds
#  (although complete target did)]

# hmmsearch --domtblout - NUCLEOTIDE VERSION!!!!!!!!!


# imports
# biopython
# os imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

#################################################################

def parse_file(HMM_std_file):
    """function to open and parse a file.
    returns """
    f = open(HMM_std_file, "r")
    # assign the file contents to the variable data
    data = f.readlines()
    # remove the \n new line and \t characters
    data1 = [line.rstrip().split("\n") for line in (data)
             if line.strip() != "" and not line.startswith("#")]
    f.close()
    return data1


def names_wanted_getter(HMM_std_file):
    """function to write the sequence not identified in the domain
    table output, but is identified as a complete macth"""
    HMM_search_data = parse_file(HMM_std_file)
    seq_names_to_get = []
    missed = ""
    for line in HMM_search_data:
        if line[0].startswith(">> "):
            missed = line[0].rstrip().replace(">> ", "")
        if "No individual domains that satisfy repor" in line[0]:
            # print("missed = ", missed)
            seq_names_to_get.append(missed)
    return seq_names_to_get


def seq_getter(filename, seq_names_to_get, out_file):
    """function to write out the seq of interest as define in the list
    from names_wanted_getter"""
    f_out = open(out_file, "w")
    for seq_record in SeqIO.parse(filename, "fasta"):
        if seq_record.id in seq_names_to_get:
            seq_record.description = ""
            SeqIO.write(seq_record, f_out, "fasta")
    f_out.close()


##########################################################################

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)

usage = """Use as follows:

Script:
# here hmmsearch looks for domian. If the whole seq matches the hmm models
# it wont put it into the tableout. So we fish for these first
# e.g. 
# >> 9_Phytophthora_irrigata_EU334634  
#   [No individual domains that satisfy reporting thresholds
# (although complete target did)]

$ python get_complete_matches_from_hmmsearch.py
-i in.fasta --HMM_std_file HMM_std.out -o outfile.fasta

requires bippython
"""

parser = OptionParser(usage=usage)

parser.add_option("-i", "--in", dest="in_file",
                  default=None,
                  help="nt_file used to generate the AA file"
                  "used for the hmmsearch")
parser.add_option("-o", "--output", dest="out_file",
                  default="nt_domains.fasta",
                  help="Output fasta filename",
                  metavar="FILE")
parser.add_option("--HMM_std_file",
                  dest="HMM_std_file",
                  default=None,
                  help="HMM_std_file --domtblout  output > HMM_std_file")


(options, args) = parser.parse_args()

filename = options.in_file
HMM_std_file = options.HMM_std_file
out_file = options.out_file

if __name__ == '__main__':
    seq_names_to_get = names_wanted_getter(HMM_std_file)
    seq_getter(filename, seq_names_to_get, out_file)






