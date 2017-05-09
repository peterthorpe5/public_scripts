#!/usr/bin/env python
# TITLE: check for illegal characters and only write out
# good sequences
# author: Peter Thorpe September 2015.
# The James Hutton Insitute, Dundee, UK.

# Imports
from optparse import OptionParser
import sys

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser

def count_letters(in_fasta, out, threshold=10):
    """this is a function check for illegal
    characters in fasta. If len is greater than threshold
    it writes out the passed seq.
    Take in a fasta, writes out a fasta"""
    threshold = int(threshold)
    aminos = set("ABRNDCQEGHIJLKMOFPSTWYVXZU")
    with open(in_fasta) as f_in:
        with open(out, 'w') as f_out:
            for title, seq in SimpleFastaParser(f_in):
                seq = seq.upper().replace("*", "")  # remove stop codons
                # TODO: Map - to X as in DIAMOND? Or remove it
                if set(seq).difference(aminos):
                    print("Check %s for %s" % (title,
                                               set(seq).difference(aminos)))
                    print(seq)
                    continue
                if len(seq)> threshold:
                    # Write out as FASTA format (no line wrapping,
                    # quicker and expect BLAST/DIAMOND not to care)
                    f_out.write(">%s\n%s\n" % (title, seq))
##    
##    f_out = open(out, 'w')
##    for seq_record in SeqIO.parse(in_fasta, "fasta"):
##        seq = str(seq_record.seq)
##        seq = seq.replace("*", "")
##        if set(seq_record.seq).difference(aminos):
##            print("Check %s for %s" % (seq_record.id,
##                                       set(seq_record.seq).difference(aminos)))
##            print(str(seq_record.seq))
##            continue
##        if len(seq_record.seq)> threshold:
##            record = SeqRecord(Seq(seq,
##            IUPAC.protein),
##            id=seq_record.id, name="",
##                description="")
##            SeqIO.write(record, f_out, "fasta")
##    f_out.close()   

if "-v" in sys.argv or "--version" in sys.argv:
    print ("v0.1.0")
    sys.exit(0)


usage = """Use as follows:
python Check..py -i nr.faa -o passed.faa -m 10

checks for illegal characters in a fasta file.
Only writes out if they are allowed.
"""

parser = OptionParser(usage=usage)
parser.add_option("-i", "--in",
                  dest="in_fasta",
                  default="nr.faa",
                  help="in filename (fasta file)",
                  metavar="FILE")

parser.add_option("-m", "--min_len",
                  dest="min_len",
                  default=10,
                  help="the min length of seq to return. " +
                  "Any seq less than this are not returned " +
                  "Default = 10")

parser.add_option("-o", "--output",
                  dest="out_file",
                  default="passed.faa",
                  help="Output filename (fasta file)",
                  metavar="FILE")

(options, args) = parser.parse_args()
outfile = options.out_file
in_fasta = options.in_fasta
min_len = options.min_len

if __name__ == '__main__':
    count_letters(in_fasta, outfile, min_len)
    print ('done')


