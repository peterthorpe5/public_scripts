#!/usr/bin/env python
# TITLE: get to cov of canu assembly for blobtools
# author: Peter Thorpe September 2015.
# The James Hutton Insitute, Dundee, UK.

# Imports
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from optparse import OptionParser


def get_cov(in_fasta, out):
    """this is a function get the cov of canu pacbio assembly.
    Take in fasta, return in file:
    name\tcov """
    f_out = open(out, 'w')
    for seq_record in SeqIO.parse(in_fasta, "fasta"):
        elements = seq_record.description.split()
        for data in elements:
            # e.g. >tig00000003 len=496087 reads=1180 covStat=487.14
            if data.startswith("covStat="):
                cov = data.split("covStat=")[1]
                out_data = "%s\t%s\n" % (seq_record.id, cov)
                f_out.write(out_data)
                continue

    f_out.close()

if "-v" in sys.argv or "--version" in sys.argv:
    print ("v0.0.1")
    sys.exit(0)


usage = """Use as follows:
python canu..py -i contigs.fasta -o contigs.cov

get the covergare of a pacbio - canu assembly
"""
parser = OptionParser(usage=usage)
parser.add_option("-i", "--in",
                  dest="in_fasta",
                  default="nr.faa",
                  help="in filename (fasta file)",
                  metavar="FILE")


parser.add_option("-o", "--output",
                  dest="out_file",
                  default="passed.faa",
                  help="Output filename (name and cov)",
                  metavar="FILE")

(options, args) = parser.parse_args()
outfile = options.out_file
in_fasta = options.in_fasta

if __name__ == '__main__':
    get_cov(in_fasta, outfile)

