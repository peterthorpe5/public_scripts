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
                 outfile, upstream, into_TSS, NNN):
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
    bind_out = outfile.split(".fa")[0] + "_%dnt_upstream_%d_into_TSS.fasta" % (upstream,
                                                                               into_TSS)
    bind_out_fa = open(bind_out, "w")
    upstream = int(upstream)
    NNN_reject_count = 0
    starting_UTR_count = 0
    fa_out_count = 0
    missing = 0
    with open(gff, "r") as f_handle:
        for line in f_handle:
            line = check_line(line)
            if not line:
                continue
            scaff, source, feature, start, stop, score, \
            direction, frame, gene_info = split_line(line)
            seq_record = genome_database[scaff]
            starting_UTR_count += 1
            if direction == "+":
                UTR = seq_record.seq[start:stop]
                bind_seq = seq_record.seq[(start - upstream):(start + into_TSS)]
                new_co_start = start - upstream
                new_co_stop = start + into_TSS
            if direction == "-":
                UTR = reverse_complement(seq_record.seq[start:stop])
                bind_seq = reverse_complement(seq_record.seq[(stop - into_TSS)
                                                             :(stop + upstream)])
                new_co_start = stop + upstream
                new_co_stop = stop - into_TSS
            coordinate_info = "\tScaffold: %s UTR_and_TSS: %d:%d  Coding_direction: %s  " % (scaff,
                                                                                 start,
                                                                                 stop,
                                                                                 direction)
            fastaextra = "returning: %d upstream of TSS and %d into UTR and or gene:" % (upstream,
                                                                                     into_TSS)
            new_coord = "  %d: %d" % (new_co_start, new_co_stop)
            description = coordinate_info + fastaextra + new_coord
            outstr = ">%s\n%s\n" % (gene_info, UTR)
            bind_str = ">%s_TSS%s\n%s\n" % (gene_info,
                                            description,
                                            bind_seq)
            if len(UTR) > min_length and len(UTR) < Max_length:
                f_out.write(outstr)
                if (NNN*"N") in bind_seq:
                    NNN_reject_count += 1
                    continue  #  we dont want NNNs
                if len(bind_seq) >= upstream: 
                    bind_out_fa.write(bind_str)
                    fa_out_count += 1
            else:
                missing += 1
    print("%d number of genes were rejected due to NNN in seq" % NNN_reject_count)
    print("we had %d UTR and TSS predictions" % starting_UTR_count)
    print("we have outputted %d fasta sequences" % fa_out_count)
    print("missing due to length problems %d" % missing)
    f_out.close()
    bind_out_fa.close()


#############################################################################
#to run it:


usage = """Use as follows:

python GFF_to_fasts.py --gff transtart.gff -m min len -x max len
        -g genome.fasta -o UTR.fasta

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
                  default=4000,
                  help="max_length of seq to return",
                  metavar="FILE")

parser.add_option("-i", "--into_TSS",
                  dest="into_TSS",
                  default=0,
                  help="into_TSS of seq to return to the binding theory",
                  metavar="FILE")

parser.add_option("-u", "--upstream",
                  dest="upstream",
                  default=20,
                  help="upstream of TSS to return",
                  metavar="FILE")

parser.add_option("-n", "--NNN",
                  dest="NNN",
                  default=300,
                  help="number of NNN in a row to be rejected from the output",
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
# - upstream
upstream = int(options.upstream)
# -x
max_length = int(options.max_length)
# --min_length
min_length = int(options.min_length)
# into_TSS
into_TSS = int(options.into_TSS)
# -n
NNN = int(options.NNN)

#######################################################################
# Run as script
if __name__ == '__main__':
    # no logging for this.
    if not os.path.isfile(genome):
        sys.exit("Input genome file not found: %s" % genome)
    if not os.path.isfile(gff):
        sys.exit("Input gff file not found: %s" % gff)
    gff_to_fasta(gff, genome, min_length,
                 max_length, outfile, upstream,
                 into_TSS, NNN)
