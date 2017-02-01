#!/usr/bin/env python
#os importsimport os
from sys import stdin,argv
import sys
from optparse import OptionParser

#Title:
#script to reformt fasta ID names.

#############################################################################
#functions
def reformat_fasta_name(filename, prefix, out_prefix="seq"):
    "this function re-write a file as a fasta file but with altered names"
    outfile = filename.split(".fa")[0]+"_alt.fasta"
    f= open(outfile, 'w')
    f_in = open(fasta, "r")
    prefix = str(prefix)
    name = []
    for seq_record in SeqIO.parse(filename, "fasta"):
        name.append(seq_record.id)
        #remove the read prefix to get to uniq names
        seq_record.id = str(seq_record.id).replace("prefix", out_prefix)
        seq_record.id= seq_record.id.replace(":", "")
        SeqIO.write(seq_record, f, "fasta")
    name_out = open("original_read_names.txt", "w")
    #file to keep track of the original names if we need them
    for i in sorted(name):
       name_out.write(i)
    name_out.close()
    f.close()
    f_in.close()




#################################################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

$ python re_format_fasta.py -f seq.fasta --prefix @M01157:20:000000000-D07KA:

script to reformt fasta names. The read name can bre

requires Biopython
"""

parser = OptionParser(usage=usage)

parser.add_option("-f","--fasta", dest="fasta", default=None,
                  help="fasta file to have names altered")

parser.add_option("-p", "--prefix", dest="prefix", default=None,
                  help="prefix to alter the id names",
                  metavar="FILE")
parser.add_option("-o", "--out_prefix", dest="out_prefix", default="seq",
                  help="prefix to the output filenames")


(options, args) = parser.parse_args()

fasta = options.fasta
prefix = options.prefix
out_prefix = options.out_prefix

#run the program

#biopython imports
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
reformat_fasta_name(fasta, prefix, out_prefix)

