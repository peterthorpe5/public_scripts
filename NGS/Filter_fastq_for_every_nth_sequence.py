
######################################################################
# Title: Script to bring back every 10000 sequences from a fastq file#
######################################################################

""" This script uses Biopython magic to iterate through a fastq file
and return every 10000 sequences.

why? : Sometimes it is good to work with a subset of data for time
and CPU reasons, as the whole data set can take much longer"""

#imports
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os
from sys import stdin,argv
import sys
from optparse import OptionParser


#####################################################################
def filter_my_fastq_file (in_fastq, number_of_seq, out_fastq):
    #open the fastq file
    number_of_seq = int(number_of_seq)
    in_file = open(in_fastq)
    #creat a new fastq file to write to
    out_file = open(out_fastq, "w")
    # enumerate is a way of counting i
    # iterate through the fatsq file
    for i, (title, seq, qual) in enumerate(FastqGeneralIterator(in_file)):
        #python magic to identify every number_of_seq "loop"
        if i % number_of_seq ==0:
            #write this to a file
            out_file.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    out_file.close()
    in_file.close()
    return True

#####################################################################



#print filter_my_fastq_file(argv[1],argv[2))


if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.2"
    sys.exit(0)


usage = """Use as follows:


python Filter....py -i infile.fastq -n every_n_number_of_seq -o outfile
"""

parser = OptionParser(usage=usage)
parser.add_option("-i", "--in", dest="in_file", default=None,
                  help="Output filename",
                  metavar="FILE")
parser.add_option("-n", "--filter_every", dest="filter", default=None,
                  help="a file with a list of gene names to get the upstream regions for")

parser.add_option("-o", "--output", dest="out_file", default=None,
                  help="Output filename",
                  metavar="FILE")





(options, args) = parser.parse_args()

infile = options.in_file
filter_n = options.filter
outfile = options.out_file

filter_my_fastq_file(infile, filter_n, outfile)
                           
