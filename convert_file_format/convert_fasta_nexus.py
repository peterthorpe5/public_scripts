# script to reformat fasta to nexus format
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
#os imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser

def convert_file(in_file, out_file):
    alignment = AlignIO.read(open(in_file), "fasta", alphabet=Gapped(IUPAC.protein))
    g = open(out_file, "w")
    g.write (alignment.format("nexus"))

##records = SeqIO.parse("G_B_M_seed_chopped002.fasta", "fasta")
###count = SeqIO.write(records, "G_B_M_seed_chopped002.txt", "phylip")
##count = SeqIO.write(records, "G_B_M_seed_chopped002.nexus", "nexus")
##print "Converted %i records" % count

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

$ convert_fasta_nexus -i in.fasta -o outfile.nexus

requires bippython
"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="in_file", default=None,
                  help="protein_file")
parser.add_option("-o", "--output", dest="out_file", default=None,
                  help="Output filename",
                  metavar="FILE")
parser.add_option("-d",  dest="DNA_file", default=None,
                  help="DNA_infile")




(options, args) = parser.parse_args()

in_file = options.in_file
out_file = options.out_file

(options, args) = parser.parse_args()

convert_file(in_file, out_file)
