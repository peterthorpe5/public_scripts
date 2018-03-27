#!/usr/bin/env python3
# title: GFF to fasta
# (c) The James Hutton Institute 2018
# Author: Peter Thorpe. The James Hutton Institute, Uk
# imports
import sys
import os
from optparse import OptionParser
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

print("""Warning: This is made for a specific purpose. We specific data
      If you are using it, check it works for you!!""")


def index_genome_file(genome):
    """index the genome file
    return a seq_record object that can be accessed
    in a dictionary like manner"""
    genome_database = SeqIO.index(genome, "fasta")
    return genome_database


def check_line(line):
    """checks if line is ok
    Starts with a comment of blank ine = not ok. """
    if line.startswith("#"):
        return False
    if not line.strip():
        return False
    return line


def split_line(line):
    """split the gff line
    Takes in a gff line. returns the elements, as str or int"""
    warning_list = ["gene", "exon", "intron"]
    assert len(line.split("\t")) == 9 , "GFF fields wrong length should be 9"
    scaff, source, feature, start, stop, score, \
            direction, frame, gene_info = line.split("\t")
    if feature in warning_list:
        print("This script is for transtart output only. Be carfeful!!!")
    gene_info = gene_info.rstrip()
    start = int(start)
    stop = int(stop)
    return scaff, source, feature, start, stop, score, \
            direction, frame, gene_info


def gff_to_fasta(gff, genome, min_length, Max_length,
                 outfile):
    """take in gff file. Gets the seq defined by the gff coords.
    If negative direction coding, the reverse complement is generated.
    A min length of seq to return and max len is applied to remove seq
    less than, for example 3 which cant be real and less that e.g., 25k
    which will be flase positives and not informative in downstream analysis
    """
    print("Indexing the genome")
    min_length = int(min_length)
    genome_database = index_genome_file(genome)
    print("Now iterating through the GFF. Assume it is sorted")
    f_out = open(outfile, "w")
    with open(gff, "r") as f_handle:
        for line in f_handle:
            line = check_line(line)
            if not line:
                continue
            scaff, source, feature, start, stop, score, \
            direction, frame, gene_info = split_line(line)
            seq_record = genome_database[scaff]
            if direction == "+":
                seq_record.seq = seq_record.seq[start:stop]
            if direction == "-":
                seq_record.seq = reverse_complement(seq_record.seq[start:stop])
            outstr = ">%s\n%s\n" % (gene_info, seq_record.seq)
            if len(seq_record.seq) > min_length:
                f_out.write(outstr)
    f_out.close()


#############################################################################
#to run it:


usage = """Use as follows:

python GFF_to_fasts.py --gff transtart.gff -g genome.fasta -o UTR.fasta

"""

parser = OptionParser(usage=usage)

parser.add_option("--gff", dest="gff",
                  default=None,
                  help="the predicted coordinate for the cds predictions .gff",
                  metavar="FILE")

parser.add_option("-g", "--genome",
                  dest="genome",
                  default=None,
                  help="the genome sequence.",
                  metavar="FILE")

parser.add_option("-m", "--min_length",
                  dest="min_length",
                  default=8,
                  help="min_length of seq to return",
                  metavar="FILE")

parser.add_option("-x", "--max_length",
                  dest="max_length",
                  default=8000,
                  help="max_length of seq to return",
                  metavar="FILE")

parser.add_option("-o", "--out", dest="outfile",
                  default="transtart.UTR.fasta",
                  help="Output filename (transtart.UTR.fasta)",
                  metavar="FILE")

(options, args) = parser.parse_args()

#-g
genome = options.genome
#--gff
gff = options.gff
#-o
outfile= options.outfile

if not os.path.isfile(genome):
    sys.exit("Input BAM file not found: %s" % genome)

#######################################################################
# Run as script
if __name__ == '__main__':
    # no logging for this.
    gff_to_fasta(gff, genome, options.min_length,
                 options.max_length, outfile)
    


