#!/usr/bin/env python
#os importsimport os
from sys import stdin,argv
import sys
from optparse import OptionParser

#Title:
#script to reformt GFF. for MCScanX

#############################################################################
#functions

def reformat_gff_column_9(gff, species, out, mcscan):
    "function to reformat name. Remove pipes and make names shorter"
    f_in = open(gff, "r")
    f_out = open(out, 'w')

    for line in f_in:
        if line.startswith("#"):
            continue
        assert len(line.split("\t")) ==9 ,"GFF fields wrong length should be 9"
        scaf,aug,cds_type,start,stop,e,f,g,gene_info = line.split("\t")
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
        #print "gene_info =", gene_info
        #print "start = ", start
        #print "stop =", stop
        if not cds_type == "gene":
            continue
        data = "%s%s\t%s\t%s\t%s\n" %(species,scaf,gene_info,\
                                    str(start), str(stop))
        f_out.write(data)
    f_out.close()
    f_in.close()


#################################################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

$ python re_format_gff_mcscan.py --gff augustus.gff -s Mc -o Gff_for_mcscan


"""

parser = OptionParser(usage=usage)

parser.add_option("--gff", dest="gff", default=None,
                  help="hintsfile",
                  metavar="FILE")
parser.add_option("-m", "--mcscan", dest="mcscan", default=False,
                  help="specific formatting for Mcscan. Default is false, change to True if you need this. ",
                  metavar="FILE")
parser.add_option("-s", "--species", dest="species", default="Sp",
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
reformat_gff_column_9(gff, species, out, mcscan)

