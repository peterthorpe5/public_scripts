#!/usr/bin/env python
# TITLE: check for illegal characters and only write out
# good sequences
# author: Peter Thorpe September 2015.
# The James Hutton Insitute, Dundee, UK.

# Imports
from optparse import OptionParser
import sys
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser

def check_illegal_character(sequence):
    """function to check if seq contains anything
    other than ATGC"""
    if set(seq).difference(nucletides):
        print("Check %s for %s" % (title,
                                  set(seq).difference(nucletides)))
        print("not adding this to database")
        illegal_character_count += 1
    

def count_letters(in_fasta, out, threshold=10):
    """this is a function check for illegal
    characters in fasta. If len is greater than threshold
    it writes out the passed seq.
    Take in a fasta, writes out a fasta"""
    threshold = int(threshold)
    nucletides = set("ATCG")
    names = set([])
    name_seq_dict = dict()
    illegal_character_count = 0
    duplciated_seq_count = 0
    f_out = open(out, 'w')
    with open(in_fasta) as f_in:
        for title, seq in SimpleFastaParser(f_in):
            if seq in name_seq_dict:
                print("duplicate seq found %s\t matches %s" % (title,
                                                               name_seq_dict[seq]))
                duplciated_seq_count += 1
                name = "_" + title.replace("Phytophthora", "P.")
                if set(seq).difference(nucletides):
                    print("Check %s for %s" % (title,
                                               set(seq).difference(nucletides)))
                    print("not adding this to database")
                    illegal_character_count += 1
                else:
                    name_seq_dict[seq] += name
            else:
                name_seq_dict[seq] = title
            if title in names:
                print("duplicate name found %s" % title)
                names.add(title)
    for keys, vals in name_seq_dict.items():
        sequence = str(keys)
        entry = ">%s\n%s\n" % (vals, sequence)
        if set(sequence).difference(nucletides):
            print("Check %s for %s" % (vals,
                                       set(sequence).difference(nucletides)))
            print("not adding this to database")
            continue
        else:
            f_out.write(entry)
    f_out.close()
    print("%d\tsequences with illegal characters" % illegal_character_count)
    print("%d\tduplciated sequences" % duplciated_seq_count)


if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.1.0")
    sys.exit(0)


usage = """Use as follows:
python Check..py -i nt.fa -o passed.fa -m 10

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

