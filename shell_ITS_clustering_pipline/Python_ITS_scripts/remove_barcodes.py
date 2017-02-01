#!/usr/bin/env python
# os imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser

# Title:
# script to trim barcodes or primer seq.

########################################################
# functions

def trim_seq(seq, barcode, left_primer, right_primer):
    """function to trim the seq at given features"""
    seq = seq.upper()
    seq = seq.replace("\n" , "")
    left_primer = int(left_primer)
    right_primer = int(right_primer)
    # remove the left and right barcode
    seq = seq[barcode:len(seq)-barcode]
    if left_primer or right_primer:
        #if given left and right primer. Slice these off
        seq = seq[left_primer:len(seq)-right_primer]
    # index the start of the ITS seq, slice there.
    # if use_ITS:
        # seq = seq[seq.index(ITS):]
        # assert seq[:6] =="CCACAC"
    return seq+"\n"

def reformat_fasta_name(fasta, barcode, left_primer, \
                        right_primer, use_ITS,\
                        ITS, out):
    """function to retun the fa file with barcodes
    and primers removed"""
    f= open(out, 'w')
    f_in = open(fasta, "r")
    ITS = str(ITS)
    left_primer = int(left_primer)
    right_primer = int(right_primer)
    barcode = int(barcode)
    seq = ""
    for line in f_in:
        if line.startswith(">"):
            if seq =="":
                # first id. write it. 
                f.write(line.split(" ")[0]+"\n")
                continue
            else:
                # we have hit a new seq record. Lets wirte
                # out details. 
                sliced_seq = trim_seq(seq, barcode, \
                                      left_primer, \
                                      right_primer)
                f.write(sliced_seq)
                # reset the variable
                seq = ""
                # write out the new name
                f.write(line.split(" ")[0]+"\n")
        else:
            seq = seq+line


###########################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

$ python re_format_fasta.py -f seq.fasta --barcode 6 --left --right

default barcode length is 6. Change by passing an option

--left is the length of the left primer  - 0 by default
--right is the length of the right primer = 0 by default

GGAAGGTGAAGTCGTAACAAGG  = Primer ITS6 fwd

start of ITS = CCACAC

the script indexes for this instead.


requires Biopython
"""

parser = OptionParser(usage=usage)

parser.add_option("-f","--fasta", dest="fasta", default=None,
                  help="fasta file to have names altered")

parser.add_option("-b", "--barcode", dest="barcode", default=6,
                  help="length of barcodes",
                  metavar="FILE")


parser.add_option("-r", "--right", dest="right_primer", default=False,
                  help="length right primer",
                  metavar="FILE")

parser.add_option("-l", "--left", dest="left_primer", default=False,
                  help="length of left primer",
                  metavar="FILE")

parser.add_option("--use_ITS", dest="use_ITS", default=False,
                  help="index and chop at a given ITS seq??",
                  metavar="FILE")

parser.add_option("--ITS", dest="ITS", default="CCACAC",
                  help="length right primer",
                  metavar="FILE")

parser.add_option("-o", "--out", dest="out", default=None,
                  help="output filenames")


(options, args) = parser.parse_args()

fasta = options.fasta
barcode = options.barcode
left_primer = options.left_primer
right_primer = options.right_primer
ITS = options.ITS
use_ITS = options.use_ITS
out = options.out

#r un the program

# print "ITS = ", ITS

# biopython imports
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
reformat_fasta_name(fasta, barcode, left_primer, right_primer,\
                    use_ITS, ITS, out)

