#!/usr/bin/env python
#author: Peter Thorpe September 2016. The James Hutton Insitute,Dundee,UK.

#Title:
#script to get the gene columns only from GFF"

#imports
import os
import sys
from sys import stdin,argv
import sys
import datetime
from optparse import OptionParser


###########################################################################


def write_out_ITS_GFF(gff, out):
    """function parse and print GFF lines that
    correspond to gene only """
    gff_file = open(gff, "r")
    out_file = open(out, "w")
    for line in gff_file:
        if line.startswith("#"):
            continue
        assert len(line.split("\t")) ==9 ,"GFF fields wrong length should be 9"
        scaffold,aug,cds_type,start,stop,e,f,g,gene_info = line.split("\t")
        if cds_type =="gene":
            out_file.write(line)
    gff_file.close()
    out_file.close()
        
       


###########################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print ("v0.0.1")
    sys.exit(0)


usage = """Use as follows:
Title:
script to get the gene columns only from GFF

$ get_genes_from_GFF.py --gff gff.out -o out.gff

"""

parser = OptionParser(usage=usage)

parser.add_option("-g", "--gff", dest="gff", default=None,
                  help="predicted gene in gff3 format",
                  metavar="FILE")

parser.add_option("-o", "--out_file", dest="out_file",
                  default="ITS_GFF.out",
                  help="output line corresponding to genes only.")


(options, args) = parser.parse_args()


gff = options.gff
out_file = options.out_file



#run the program

if not os.path.isfile(gff):
    print("sorry, couldn't open the file: " + ex.strerror + "\n")
    print ("current working directory is :", os.getcwd() + "\n")
    print ("files are :", [f for f in os.listdir('.')])
    sys_exit("\n\nInput blast file not found: %s" % gff)

# call the top function    
write_out_ITS_GFF(gff, out_file)


