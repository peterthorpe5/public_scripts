#!/usr/bin/env python
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
from sys import stdin,argv
import sys
from optparse import OptionParser

def seq_getter(filename_in, wantedfile, threshold, outfile):
    "script to gt sequences of intereste from a file of wanted genes"

    f= open(outfile, 'w')
    wanted = open(wantedfile, "r")

    names = wanted.readlines()
    #print names
    wanted_data = [line.replace("\t", "").rstrip("\n") for line in names
              if line.strip() != ""]
    name_set = set([])
    for i in wanted_data:
        if not i.startswith("#"):
            i = i.rstrip()
            name_set.add(i)
    #print wanted_data

    cds_database = SeqIO.index(filename_in, "fasta")
    #record = SeqIO.read(filename, "fasta")
    for i in name_set:
        if "\r\n" in i:
            i = i.replace("\r\n","")
        #print i
        seq_record = cds_database[i]
        dashes = seq_record.seq.count("-")
        print 100*(float(dashes)/len(seq_record.seq))
        if 100*(float(dashes)/len(seq_record.seq)) < int(threshold):
            #print 'boomshanka'
            SeqIO.write(seq_record, f, "fasta")
    f.close()
    return True

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

$ python get_seq_filter_by_dashes.py -i in.fasta -t threshold_percentage_of_dahses_per_seq -o outfile.fasta

if you wish you can pass it a file of wanted gene names to further refine those which you are working with

requires bippython
"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="in_file", default=None,
                  help="domain file")
parser.add_option("-o", "--output", dest="out_file", default=None,
                  help="Output filename",
                  metavar="FILE")
parser.add_option("--wanted_genes", dest="wanted_gene", default=None,
                  help="file with a list of wanted_gene names",
                  metavar="FILE")

parser.add_option("-t",  dest="threshold", default=None,
                  help="threshold_percentage_of_dahses_per_seq")




(options, args) = parser.parse_args()

in_file = options.in_file
threshold = options.threshold
wanted_gene = options.wanted_gene
out_file = options.out_file

(options, args) = parser.parse_args()



seq_getter(in_file, wantedfile, threshold, out_file)

#seq_getter('assembly2_scaffolds.fasta',\
           #'scaffold318.fasta')
print 'done'

