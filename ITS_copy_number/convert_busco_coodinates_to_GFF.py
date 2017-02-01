#!/usr/bin/env python
#author: Peter Thorpe September 2016. The James Hutton Insitute,Dundee,UK.

#Title:
#script to convert the BUSCO coordinates to GFF for bedtools

#imports
import os
import sys
from sys import stdin,argv
import sys
import datetime
from optparse import OptionParser


###########################################################################

def parse_tab_outfile(busco):
    """read in the busco tab file. Reads whole file into memory.
    returns a list, one list item per busco hit.
    """
    with open(busco) as file:
        return file.read().split("\n")

def write_out_ITS_GFF(busco, prefix, out): # this is a long function
    """function to write out the busco hits in a GFF3
    like manner. """

    try:
        busco_hits = parse_tab_outfile(busco)
    except:
        raise ValueError("something wrong with gff in file")
    GFF_out = open(out, "w")
    
    ############################################################################
    EOG_counter = 0
    for i in busco_hits:
        if i.startswith("#"): #allows line to have comment.
            continue
        if not i.strip():
            continue #if the last line is blank
        EOG_counter = EOG_counter+1
        EOG, scaff, start, stop = i.split("\t")
        # check start is less than stop. Fix if not. 
        if int(start) > int(stop):
            temp_start = stop
            temp_stop = start
            start = temp_start
            stop = temp_stop
        if int(start) == 0:
            start = "1"
        line = "%s\t%s_Busco_gene\t%d\t%s\t%s\t.\t+\t.\t%s\n" %(scaff,\
                        prefix,EOG_counter, start, stop, EOG)
        GFF_out.write(line)
        
        

    #close the write file
    GFF_out.close()         


###########################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print ("v0.0.1")
    sys.exit(0)


usage = """Use as follows:
Title:
script to convert the BUSCO coordinates to GFF for bedtools

$ convert_busco_coodinates_to_GFF.py -b ./run_busco/coordinates_busco --prefix info_line -o out.gff

Note:
coloumns 6,7 and 8 are 'made' up for the purpose of this task.
busco hits on the negative strand will be inverted so the
start is always less than the end coordinate.

"""

parser = OptionParser(usage=usage)

parser.add_option("-b", "--busco", dest="busco", default=None,
                  help="the tab out file from the busco search",
                  metavar="FILE")
parser.add_option("--prefix", dest="prefix",
                  default="temp_name",
                  help="name for column 2 in GFF. Best to "
                  " use the origin of the data")

parser.add_option("-o", "--out_file", dest="out_file",
                  default="ITS_GFF.out",
                  help="outfile for the busco regions in GFF format")


(options, args) = parser.parse_args()


busco = options.busco
prefix = options.prefix
out_file = options.out_file



#run the program

if not os.path.isfile(busco):
    print("sorry, couldn't open the file: " + ex.strerror + "\n")
    print ("current working directory is :", os.getcwd() + "\n")
    print ("files are :", [f for f in os.listdir('.')])
    sys_exit("\n\nInput busco file not found: %s" % gff)

# call the top function    
write_out_ITS_GFF(busco, prefix, out_file)


