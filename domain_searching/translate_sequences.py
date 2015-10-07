#!/usr/bin/env python

#Title: script to translate nt sequences in to Amino acids.

#imports
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
import os
from sys import stdin,argv
import sys
from optparse import OptionParser

###############################################################################################################

def name_set(in_names):
    "function to return a set of names, to remove duplicates etc.."
    name_set = set([])
    for i in in_names:
        if not i.startswith("#"):
            name_set.add(i)
    return name_set


def translate_fasta(in_file, length, wanted, not_wanted, out):
    """ function to return the translated sequence given a length threshold and
    possible lists of specifically wanted or unwanted gene names"""
    f= open (out,'w')
    count = 0
    #get the wanted gene name list
    wanted_name_set = set([])
    if wanted:
        names = wanted.readlines()
        wanted_data = [line.replace("\t", "").rstrip("\n") for line in names
                  if line.strip() != ""]
        #call function
        wanted_name_set = name_set(wanted_data)
    not_wanted_name_set = set([])    
    #get the unwanted gene name list
    if not_wanted:
        names = not_wanted.readlines()
        not_wanted_data = [line.replace("\t", "").rstrip("\n") for line in names
                  if line.strip() != ""]
        #call function
        not_wanted_name_set = name_set(not_wanted_data)
        
    for seq_record in SeqIO.parse(in_file, "fasta"):
        #print seq_record.id
        if len(seq_record.seq)> int(length):
            if wanted:
                if seq_record.id in wanted_name_set:
                    seq_record.seq = eq_record.seq.translate()
                    SeqIO.write(seq_record, f, "fasta")
            if not_wanted:
                if not seq_record.id in not_wanted_name_set:
                    seq_record.seq = eq_record.seq.translate()
                    SeqIO.write(seq_record, f, "fasta")
            else:
                seq_record.seq = seq_record.seq.translate()
                SeqIO.write(seq_record, f, "fasta")  
        count +=1
    print ("I have translated %d sequences" % (count))
    return True

###############################################################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

Translates nt to AA sequences

$ python translate_sequences.py -i in.fasta -l min_length_of_seq (default(3)) --not_wanted --wanted -o out.fasta

script either reformats badly formated fasta file. Within reason. Ie. Word documents will still break it.

if lenght of seq is longer than -l default 3 - writes to file.

can filter fasta by giving it a list of --not_wanted or --wanted names.

"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="in_file", default="all_nt.fasta",
                  help="current fasta you want to reformat")

parser.add_option("--wanted", dest="wanted", default=False,
                  help="a text file with a list of wanted gene names found in the fasta file",
                  metavar="FILE")

parser.add_option("--not_wanted", dest="not_wanted", default=False,
                  help="a text file with a list of not_wanted gene names found in the fasta file",
                  metavar="FILE")

parser.add_option("-l", "--lenth", dest="length", default="3",
                  help="Output filename",
                  metavar="FILE")
parser.add_option("-o", "--out", dest="out", default=None,
                  help="Output filename",
                  metavar="FILE")




(options, args) = parser.parse_args()

in_file = options.in_file
not_wanted = options.not_wanted
wanted = options.wanted
out = options.out
length = options.length



translate_fasta(in_file,length,wanted,not_wanted, out)
print 'done'

