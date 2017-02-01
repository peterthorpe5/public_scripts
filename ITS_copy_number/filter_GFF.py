#!/usr/bin/env python
#author: Peter Thorpe September 2016. The James Hutton Insitute,Dundee,UK.

#Title:
#script to reduce redundancy in GFF file)"

#imports
import os
import sys
from sys import stdin,argv
import sys
import datetime
from optparse import OptionParser


###########################################################################

def parse_tab_outfile(blast):
    """read in the blast tab file. Reads whole file into memory.
    returns a list, one list item per blast hit.
    """
    with open(blast) as file:
        return file.read().split("\n")

def split_gff_line(line):
    """function to split the gff line into its components"""
    assert len(line.split("\t")) ==9 ,"GFF fields wrong length should be 9"
    scaf, genome, hit_number, start, stop, dot, direction, \
          dot2, description = line.split("\t")
    return scaf, genome, hit_number, start, stop, dot, direction, \
          dot2, description 
    


def write_out_ITS_GFF(gff, out): # this is a long function
    """function to write out the ITS blast hits in a GFF3
    like manner. """
    #this is out list of so called top/ longest matches which we will
    #append/remove as applicable
    merged_blast_hits = []
    blast_hit_to_info_dict = dict()

    best_stop = 0
    best_start = 0
    last_scaffold = "tmp"
    last_hit_number = ""
    
    # call function to get list of blast hits.
    try:
        blast_hits = parse_tab_outfile(gff)
    except:
        raise ValueError("something wrong with gff in file")
    GFF_out = open(out, "w")
    
    ############################################################################
    for i in blast_hits:
        #print "best_start, best_stop", best_start, best_stop
        if i.startswith("#"): #allows line to have comment.
            continue
        if not i.strip():
            continue #if the last line is blank
        
        assert len(i.split("\t"))== 9, ("""custom BLAST output?
                        not enough coloumns in gff file.""")
        #split up the gff line
        scaf, genome, hit_number, start, stop, dot, direction, \
                  dot2, description = split_gff_line(i)
        #populate the dictionaly
        blast_hit_to_info_dict[hit_number] = i

        #print "------------------"
        #print "scaf", scaf, "last_scaffold", last_scaffold
        #print "current :", i
        #print "start :", start, "stop", stop
        
        # this is the first iteration. Populate the variables. 
        if last_scaffold == "tmp":
            best_stop = int(stop)
            best_start = int(start)
            last_scaffold = scaf
            last_hit_number = hit_number
            merged_blast_hits.append(i)
            continue

        if scaf == last_scaffold:
            assert int(start) >= int(best_start), ("""your
                    gff file has not been sorted by linux sort
                    please run this command
                    cat ${genome_prefix}.ITS.GFF | sort -k1,1 -k4n -k5n
                    > sorted.gff .. error in line %s""" %i)

            #print "best_start ;", best_start, "best_stop", best_stop, "\n"

            same_hit = False
            # same scaffold. Is the hit in the same region.
            #this means it could the same blast hit region
            if best_start <= int(start) <= best_stop:
                same_hit = True
                # Does this fall within the current merged hit,
                # or does it extend the current merged hit?
                if int(stop) <= best_stop:
                    # Falls within the current merged hit, boring
                    #print ("i am ignoring: %s" %(hit_number))
                    del blast_hit_to_info_dict[hit_number]
                    # leave merged_blast_hits as it is.
                else:
                    #print("Extends the current merged hit which was %i %i" % (best_start, best_stop))
                    #print "stop" , stop, "best_stop", best_stop
                    #update the stop value
                    best_stop = int(stop)
                    #remove this current dictionary entry
 
                    # get the current merged values, want to use old name
                    scaf, genome, hit_number, start, stop, dot, direction, \
                          dot2, description = split_gff_line(merged_blast_hits[-1])
                    if not hit_number.endswith("_merged"):
                        hit_number += "_merged"
                    #use this old vlaues with the new end coordinate
                    updated_values = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(scaf,\
                                     genome, hit_number, \
                                     best_start, best_stop, \
                                     dot, direction, dot2, description)
                    
                    #update this blast enery in the 
                    blast_hit_to_info_dict[last_hit_number] = updated_values
                    #print("Merged hit now %i %i" % (best_start, best_stop))
                    # Replace last value with extended hit
                    merged_blast_hits[-1] = updated_values
            else:
                # This does not overlap, it is the first hit in a new
                # merged hit.
                #print("This is a new region (did not overlap with %i %i)" % (best_start, best_stop))
                best_stop = int(stop)
                best_start = int(start)
                merged_blast_hits.append(i)
            # Carry on to next hit on this scaffold...
            assert last_scaffold == scaf
        else:
            # new scafold
            # (results from last scaffold already saved in blast_hit_to_info_dict)
            best_stop = int(stop)
            best_start = int(start)
            last_scaffold = scaf
            last_hit_number = hit_number
            merged_blast_hits.append(i)
    print_list = []        
    for key, val in blast_hit_to_info_dict.items():
        print_list.append(val)
    #print "\n\nResults:"
    for concensus_hit in merged_blast_hits:
        GFF_out.write(concensus_hit+"\n")

    #close the write file
    GFF_out.close()         


###########################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print ("v0.0.1")
    sys.exit(0)


usage = """Use as follows:
Title:
script to reduce blast hits to consensus blast hit. This is required for
acurate genomic read coverage later on.

The BLAST should already have been perfomed:
 blastn -query ITS.fasta -db genome.fasta -outfmt 6 -out ITS_vs_geome.out
 a gff shoukld alread have been produced using 'generate_ITS_GFF.py'


$ filter_GFF.py --gff gff.out -o out.gff

Note:
coloumns 6,7 and 8 are 'made' up for the purpose of this task.
BLAST hits on the negative strand will be inverted so the
start is always less than the end coordinate.

"""

parser = OptionParser(usage=usage)

parser.add_option("-g", "--gff", dest="gff", default="outfmt6.out",
                  help="the tab out file from the BLAST search",
                  metavar="FILE")

parser.add_option("-o", "--out_file", dest="out_file",
                  default="ITS_GFF.out",
                  help="outfile for the ITS regions in GFF format")


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


