
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
def filter_my_fastq_file (in_fastq, read_names, out_fastq):
    #open the fastq file
    in_file = open(in_fastq)
    #creat a new fastq file to write to
    out_file = open(out_fastq, "w")

    wanted = open(read_names, "r")
    print ("im loading the read names")

    names = wanted.readlines()
    #print names
    wanted_data = [line.rstrip() for line in names
              if line.strip() != ""]
    name_set = set([])
    for i in wanted_data:
        if not i.startswith("#"):
            name_set.add(i)
    
    #print wanted_data
    names = set([])
    print ("I have loaded the read names")
    
    # enumerate is a way of counting i
    # iterate through the fatsq file
    for i, (title, seq, qual) in enumerate(FastqGeneralIterator(in_file)):
        if title in names:
            #write this to a file
            out_file.write("@%s\n%s\n+\n%s\n" % (title, seq.upper(), qual))
    out_file.close()
    in_file.close()
    return True

#####################################################################



#print filter_my_fastq_file(argv[1],argv[2))


if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.2"
    sys.exit(0)


usage = """Use as follows:


python get_fastq....py -i infile.fastq -n names_of_reads -o outfile
"""

parser = OptionParser(usage=usage)
parser.add_option("-i", "--in", dest="in_file", default=None,
                  help="Output filename",
                  metavar="FILE")
parser.add_option("-n", "--names", dest="filter", default=None,
                  help="a file of read name. Has to match exactly the names in fq file")

parser.add_option("-o", "--output", dest="out_file", default=None,
                  help="Output filename",
                  metavar="FILE")





(options, args) = parser.parse_args()

infile = options.in_file
filter_n = options.filter
outfile = options.out_file

filter_my_fastq_file(infile, filter_n, outfile)
                           
