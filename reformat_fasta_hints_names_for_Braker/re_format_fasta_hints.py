#!/usr/bin/env python
#os importsimport os
from sys import stdin,argv
import sys
from optparse import OptionParser

#Title:
#script to reformt scaffold names for Braker.
#remove pipes and reduce the length of the names

#############################################################################
#functions
def reformat_fasta_scaffold_name(filename, prefix, out_prefix):
    "this function re-write a file as a fasta file but with altered names"
    outfile = out_prefix+"_alt.fasta"
    f= open(outfile, 'w')
    f_in = open(fasta, "r")
    prefix = str(prefix)
    for seq_record in SeqIO.parse(filename, "fasta"):
        seq_record.id = str(seq_record.id).split("|")[0]
        seq_record.id= seq_record.id.replace("scaffold", prefix)
        seq_record.description = ""
        SeqIO.write(seq_record, f, "fasta")                    
    f.close()
    f_in.close()
    return True

def reformat_hints_scaffold_name(hints, prefix, out_prefix):
    "function to reformat name. Remove pipes and make names shorter"
    f_in = open(hints, "r")
    outfile = out_prefix+"_alt.hints"
    f_out = open(outfile, 'w')
    prefix = str(prefix)
    for line in f_in:
        if "scaffol" in line:
            #line = line.replace("|", "_")
            scaf,a,b,c,d,e,f,g,h = line.split("\t")
            scaf = scaf.split("|")[0]
            scaf = scaf.replace("scaffold", prefix)
            data = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(scaf,a,b,c,d,e,f,g,h.rstrip())

        print >>f_out,data
    f_out.close()
    f_in.close()


#################################################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

$ python re_format_fasta_hints -f genome.fasta --hints hints_file --prefix Mc

script to reformt scaffold names for Braker. It doesnt seems to like '|' or long names

requires Biopython!! python 2.7 - if you have python 3 alter this line: print >>f_out,data and add a "\n" to data
"""

parser = OptionParser(usage=usage)

parser.add_option("-f","--genome", dest="fasta", default=None,
                  help="fasta genome file to have names altered")
parser.add_option("--hints", dest="hints", default=None,
                  help="hintsfile",
                  metavar="FILE")
parser.add_option("-p", "--prefix", dest="prefix", default="Sc",
                  help="prefix to alter the full 'scaffold123|length1234' to 'prefix123'",
                  metavar="FILE")
parser.add_option("-o", "--out_prefix", dest="out_prefix", default=None,
                  help="prefix to the output filenames e.g. Out_prefix.fasta out_prefix.hints")


(options, args) = parser.parse_args()

fasta = options.fasta
hints = options.hints
prefix = options.prefix
out_prefix = options.out_prefix

#run the program

reformat_hints_scaffold_name(hints, prefix, out_prefix)
#biopython imports required for this function
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
reformat_fasta_scaffold_name(fasta, prefix, out_prefix)

