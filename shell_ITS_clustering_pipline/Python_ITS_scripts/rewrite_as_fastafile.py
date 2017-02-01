#!/usr/bin/env python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
#os imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser

def name_set(in_names):
    name_set = set([])
    for i in in_names:
        if not i.startswith("#"):
            name_set.add(i)
    return name_set
    

def reformat_as_fasta(filename, min_length_of_seq, outfile):
    "this function re-write a file as a fasta file"
    f= open(outfile, 'w')

    for seq_record in SeqIO.parse(filename, "fasta"):
        if len(seq_record.seq) >= int(min_length_of_seq):
            SeqIO.write(seq_record, f, "fasta")                    
    
    f.close()
    return True



if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

converts

$ python rewrite_as_fasta.py -i in.fasta -l min_length_of_seq (default(3)) --not_wanted --wanted -o out.fasta

script either reformats badly formated fasta file. Within reason. Ie. Word documents will still break it.



"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="in_file", default=None,
                  help="current fasta you want to reformat")
parser.add_option("-m", dest="min_length_of_seq", default=130,
                  help="min_length_of_seq")
parser.add_option("-o", "--out", dest="out", default=None,
                  help="Output filename",
                  metavar="FILE")




(options, args) = parser.parse_args()

in_file = options.in_file
min_length_of_seq = options.min_length_of_seq
out = options.out



reformat_as_fasta(in_file, min_length_of_seq, out)
print 'done'

