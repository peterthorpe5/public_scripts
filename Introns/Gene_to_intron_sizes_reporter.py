#!/usr/bin/env python
# author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

# imports
import os
import sys
from sys import stdin,argv
import sys
from optparse import OptionParser
import datetime

# Title:
# How many exons does a gene have?

#############################################################################

def split_gff_gene_names(gene_info):
    """function to return just the gene name"""
    #print gene_info
    try:
        gene = gene_info.split("ID=")[1]
    except:
        gene = gene_info
    try:
        gene = gene_info.split("Parent=")[1]
    except:
        gene = gene
    gene = gene.split(".gene")[0]
    gene = gene.split(".exon")[0]
    gene = gene.split(";")[0]
    gene = gene.split(";")[0]
    gene = gene.split(".CDS")[0]
    gene = gene.split(".t")[0]
    # some data set require this to be uncommented
    # gene = gene.split(".t")[0]
    gene = gene.rstrip("\n")
    #gene = gene.replace(".t1", "")
    return gene


def parse_gff(gff, out_file):
    """function to parse GFF and produce a scaffold_to_gene_dict
    3 dictionaries are returned:
    scaffold_to_gene_dict, gene_to_exon_count, gene_start_stop_dict"""
    f_in = open(gff, "r")
    # data to output of function
    scaffold_to_gene_dict = dict()
    gene_to_exon_count = dict()
    gene_start_stop_dict = dict()
    count = 1
    intron_size_out = open(out_file, "w")
    title_intronsizes = "#gene\tintron_number\tstart\tstop\tintron_length\n"
    intron_size_out.write(title_intronsizes)
    # iterate through gff
    for line in f_in:
        if line.startswith("#"):
            continue
        if not line.strip():
                continue  #  if the last line is blank
        scaffold, aug, cds_type, start, stop, e, f, \
                  g,gene_info = line.split("\t")
        gene = split_gff_gene_names(gene_info)
        if not gene in gene_to_exon_count:
            gene_to_exon_count[gene] = 0
        # scaffold_to_gene_dict
        if line.split("\t")[2] == "gene":
            if not scaffold in scaffold_to_gene_dict:
                scaffold_to_gene_dict[scaffold]=[gene]
            else:
                scaffold_to_gene_dict[scaffold].append(gene)
            start_stop_formatted = "%s\t%d\t%d" %(scaffold, int(start),
                                                  int(stop))
            gene_start_stop_dict[gene] = start_stop_formatted

        # gene_to_exon_count
        if line.split("\t")[2] == "intron":
            gene_to_exon_count[gene] += 1
            intron_len = int(stop) - int(start)
            outfmt = "%s\tintron_%d\t%s\t%s\t%d\n" % (gene,
                                                      gene_to_exon_count[gene],
                                                      start,
                                                      stop,
                                                      intron_len)
            intron_size_out.write(outfmt)
    # print scaffold_to_gene_dict
    f_in.close()
    print("if this crashes check exons are in your gff")
    print("if not change lines 96 'exon' to 'intron'")
    intron_size_out.close()


#################################################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:



$ python Gene_to_exon_number_reporter.py --gff ../augustus.gff3 -o out.file


"""

parser = OptionParser(usage=usage)

parser.add_option("--gff", "--GFF",
                  dest="gff",
                  default="test.gff",
                  help="gff file for predicted genes. ",
                  metavar="FILE")

parser.add_option("-o", "--out_file",
                  dest="out_file",
                  default="test_output.txt",
                  help="outfile to gene to exon counts")

(options, args) = parser.parse_args()


gff = options.gff
out_file = options.out_file



#run the program
# Run as script
if __name__ == '__main__':
    # call the main function
    if not os.path.isfile(gff):
        print("sorry, couldn't open the file: " % gff)
    parse_gff(gff, out_file)

