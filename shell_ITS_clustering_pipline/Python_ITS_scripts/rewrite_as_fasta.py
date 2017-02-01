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
    

def reformat_as_fasta(filename,length,wanted,not_wanted, outfile):
    "this function re-write a file as a fasta file"
    f= open(outfile, 'w')
    database =open("name_to_species_database.txt", "w")

        
    #print wanted_data
    allowed_list = ["A", "T", "G", "C"]
    for seq_record in SeqIO.parse(filename, "fasta"):
        seq_record.seq = str(seq_record.seq)
        des = seq_record.description
        seq_record.seq = seq_record.seq.replace("N", "")
        seq_record.seq = seq_record.seq.replace("Y", "")
        seq_record.seq = seq_record.seq.replace("S", "")
        seq_record.seq = seq_record.seq.replace("K", "")
        seq_record.seq = seq_record.seq.replace("R", "")
        seq_record.seq = seq_record.seq.replace("W", "")
        seq = seq_record.seq.replace("M", "")
        
        name_ITS = des.replace(" 	(", "_")
        name_ITS = name_ITS.replace(")", "")
        print >> f ,">%s\n%s" %(seq_record.id, seq)
        print >>database, seq_record.description
        #SeqIO.write(seq_record, f, "fasta")                    
    
    f.close()
    return True



if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

converts

$ python rewrite_as_fasta.py -i in.fasta -l min_length_of_seq (default(3)) --not_wanted --wanted -o out.fasta

script either reformats badly formated fasta file. Within reason. Ie. Word documents will still break it.

if lenght of seq is longer than -l default 3 - writes to file.

can filter fasta by giving it a list of --not_wanted or --wanted names.

"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="in_file", default=None,
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


reformat_as_fasta(in_file,length,wanted,not_wanted, out)
print 'done'

