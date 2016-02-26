#!/usr/bin/env python
#os importsimport os
from sys import stdin,argv
import sys
from optparse import OptionParser

#Title:
#script to reformt GFF. This brakes tools = pri=4;src=E "664.795;0.1;0:1"

#############################################################################
#functions

def reformat_gff_column_9(gff, split_at, out_prefix):
    "function to reformat name. Remove pipes and make names shorter"
    f_in = open(gff, "r")
    outfile = out_prefix+".gff"
    f_out = open(outfile, 'w')
    split_at = str(split_at)
    for line in f_in:
        if line.startswith("#"):
            continue
        scaf,a,b,c,d,e,f,g,h = line.split("\t")
        new_coloumn_9 = h.split(split_at)[0]
        data = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(scaf,a,b,c,d,e,f,g,new_coloumn_9.rstrip())

        print >>f_out,data
    f_out.close()
    f_in.close()


#################################################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

$ python re_format_gff.py --gff augustus.gff -s src=E (default) -o agustus_reformatted

script to reformt gff coloumn 9 as the ; formatting brakes some tools. 

"""

parser = OptionParser(usage=usage)

parser.add_option("--gff", dest="gff", default=None,
                  help="hintsfile",
                  metavar="FILE")
parser.add_option("-s", "--split_at", dest="split_at", default="src=E ",
                  help="split_at   src=E  will split this line up pri=4;src=E '664.795;0.1;0:1'",
                  metavar="FILE")
parser.add_option("-o", "--out_prefix", dest="out_prefix", default=None,
                  help="prefix to the output filenames")


(options, args) = parser.parse_args()


gff = options.gff
split_at = options.split_at
out_prefix = options.out_prefix

#run the program
reformat_gff_column_9(gff, split_at, out_prefix)

