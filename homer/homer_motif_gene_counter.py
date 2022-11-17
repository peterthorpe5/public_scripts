# Title: script to count transposon types
# author: Peter Thorpe September 2015. The James Hutton Insitute, Dundee, UK.

# imports
import random
import os
import sys
from optparse import OptionParser
import datetime
from collections import defaultdict


def parse_goi(inset, infile):
    """parse the gio file"""
    f_in = open(bed, "r")
    # iterate through 
    for line in f_in:
        if line.startswith("#"):
            continue
        if not line.strip():
                continue  #  if the last line is blank
        data = line.split()
        gene = data[0].rstrip()
        inset.add(gene)
    return inset
        
    

def parse_bed(bed):
    """function to parse the elements of the bed file
    counts the occurance of the genes to see if there are
    multiple motifs found per gene. 
    INPUT e.g.
A1CF	187	387	foxp2_motif	11.066319	-
ABL2	208	408	foxp2_motif	11.066319	+
ACTA2	618	818	foxp2_motif	11.066319	+
AFF3	89	289	foxp2_motif	11.066319	+
AKAP9	676	876	foxp2_motif	11.066319	+
ALKBH8	48	248	foxp2_motif	11.066319	-
    """
    f_in = open(bed, "r")
    # data to output of function
    gene_dict = defaultdict(int)
    gene_dict_coordinates = defaultdict(list)
    count = 1

    # iterate through bed
    for line in f_in:
        if line.startswith("#"):
            continue
        if not line.strip():
                continue  #  if the last line is blank
        gene, start, stop, motif, prob, direction = line.split("\t")
        gene_dict[gene] += 1
        start_stop = "%s:%s" % (start, stop)
        gene_dict_coordinates[gene].append(start_stop)

    return gene_dict, gene_dict_coordinates

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.3")
    sys.exit(0)


usage = """Use as follows:

$ python ...counter.py --bed result.bed -o outfile.out

"""

parser = OptionParser(usage=usage)

parser.add_option("--bed",
                  dest="bed",
                  default="example.bed",
                  help="bed file output from homer",
                  metavar="FILE")

parser.add_option("--goi",
                  dest="goi",
                  default="Genes227.txt",
                  help="GOI list",
                  metavar="FILE")

parser.add_option("-o", "--out",
                  dest="out_file",
                  default="example.counts",
                  help="outfile for the results",
                  metavar="FILE")


(options, args) = parser.parse_args()

#-o
out_file = options.out_file
#-t
bed = options.bed
#--goi
goi = options.goi

#########################################################
# run the program
# Run as script
if __name__ == '__main__':
    # call the main function
    if not os.path.isfile(bed):
        print("sorry, couldn't open the file: ", bed)
    DE = set([])
    goi_227 = set([])
    DE = parse_goi(DE, "GFP_striatum_vs_KI_striatum.GLM.edgeR.DE.LOG1.2_FDR0.05_KI_striatum_UP.subset")
    goi_227 = parse_goi(goi_227, "Genes227.txt")
    print(DE)
    print(goi_227)
    gene_dict, gene_dict_coordinates = parse_bed(bed)
    outfile = open(out_file, "w")

    title = "#gene\tFOXP2_motif_Count\tin_227_genes\tDE_up_striatum\tmotif_coordinates\n"
    outfile.write(title)
    for gene, count in gene_dict.items():
        GOI227 = "-"
        DE_found = "-"
        if gene in goi_227:
            GOI227 = "YES"
        if gene not in goi_227:
            GOI227 = "-"
        if gene in DE:
            DE_found = "YES"
        if gene not in DE:
            DE_found = "-"
        coordinates = gene_dict_coordinates[gene]
            
        outfmt = "%s\t%d\t%s\t%s\t%s\n" % (gene, count, GOI227, DE_found, str(coordinates))
        #print(outfmt)
        outfile.write(outfmt)
    outfile.close()
    
