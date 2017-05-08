#!/usr/bin/env python
# Title:
# script to reformt GFF. for MCScanX, or
# for genomic upstream regions python script
# author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

import os
from sys import stdin,argv
import sys
from optparse import OptionParser

# function

def reformat_gff_column_9(gff, species, out, mcscan):
    """function to reformat name. Remove pipes and make names shorter"""
    f_in = open(gff, "r")
    f_out = open(out, 'w')

    for line in f_in:
        if line.startswith("#"):
            continue
        assert len(line.split("\t")) ==9 ,"GFF fields wrong length should be 9"
        scaf, source, feature, start, stop, score, \
              direction, frame, gene_info = line.split("\t")
        try:
            if mcscan == True:
                scaf = scaf.split("scaffold_")[1]
        except:
            ValueError
            if mcscan == True:
                scaf = scaf.replace(species, "")
        gene_info = gene_info.replace("ID=", "").split()[0]
        gene_info = gene_info.split(".t")[0]
        gene_info = gene_info.replace(";", "")
        gene_info = gene_info.split("Note=")[0]
        
        if not feature == "gene":
            continue
        if mcscan == True:
            data = "%s%s\t%s\t%s\t%s\n" %(species,
                                          scaf,
                                          gene_info,
                                          str(start),
                                          str(stop))
        else:
            # we want to add the direction information
            data = "%s%s\t%s\t%s\t%s\t%s\n" %(species,
                                              scaf,
                                              str(start),
                                              str(stop),
                                              direction,
                                              gene_info)
        f_out.write(data)
    f_out.close()
    f_in.close()


#################################################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.2"
    sys.exit(0)


usage = """Use as follows:

if reformatting gff for mcscanx type --mcscan.
$ python re_format_gff_mcscan.py --gff augustus.gff -s Mc -o Gff_for_mcscan --mcscan

else,
this will pull the gene regions out of a gff
"""

parser = OptionParser(usage=usage)

parser.add_option("--gff", dest="gff", default=None,
                  help="hintsfile",
                  metavar="FILE")
parser.add_option("-m", "--mcscan", dest="mcscan",
                  default=False,
                  help="specific formatting for Mcscan. Default is " +
                  "false, change to True if you need this. ",
                  metavar="FILE")
parser.add_option("-s", "--species", dest="species",
                  default="",
                  help="species prefix to add into column 1 two letters",
                  metavar="FILE")
parser.add_option("-o", "--out", dest="out", default=None,
                  help="output filenames")


(options, args) = parser.parse_args()


gff = options.gff
mcscan = options.mcscan
species = options.species
out = options.out

#run the program
# Run as script
if __name__ == '__main__':
    if not os.path.isfile(gff):
        print("sorry cannot find you %s file" % gff)
        os._exit(0)
    print ("warning: GFF files vary! make sure you check the file before using")
    reformat_gff_column_9(gff, species, out, mcscan)

