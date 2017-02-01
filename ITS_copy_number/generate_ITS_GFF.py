#!/usr/bin/env python
#author: Peter Thorpe September 2016. The James Hutton Insitute,Dundee,UK.

#Title:
#script to generate gff for ITS region BLAST hits)"
# The BLAST should already have been perfomed:
# blastn -query ITS.fasta -db genome.fasta -outfmt 6 -out ITS_vs_geome.out

#imports
import os
import sys
from sys import stdin,argv
import sys
import datetime
from optparse import OptionParser


###########################################################################

def parse_blast_tab_outfile(blast):
    """read in the blast tab file. Reads whole file into memeroy.
    returns a list, one list item per blast hit.
    """
    with open(blast) as file:
        return file.read().split("\n")


def get_unique_hits(temp_blast_hits):
    """function to remove duplicate hits"""
    identified_data_set = set([])
    blast_hits = []
    for result in temp_blast_hits:
        blast_line = result.split("\t")
        if len(blast_line) < 8:
            continue
        scaffold = blast_line[1]
        start = blast_line[8]
        stop = blast_line[9]
        hit = "%s\t%s\t%s" %(scaffold, start, stop)
        hit_rev = "%s\t%s\t%s" %(scaffold, stop, start)
        #print hit
        # chexck to see if this start stop scaff location
        # has already been found.
        if hit not in identified_data_set:
            identified_data_set.add(hit)
            identified_data_set.add(hit_rev)
            blast_hits.append(result)
    #print "temp blast hits", blast_hits
    best_blast_hits = get_representative_blast_hit(blast_hits)
    #print "best_blast_hits :", best_blast_hits
    return blast_hits
    
    
def spit_blast_data(i, blast_count):
    """function to split up the blast hits
    and return \t formatted data. checks start < stop
    , alters these if need be..."""
    # split the blast line and assign the feilds respectively   
    queryId, subjectId, percIdentity, alnLength,mismatchCount,\
             gapOpenCount, queryStart, queryEnd, subjectStart, \
             subjectEnd, eVal, bitScore = i.split("\t")
    #reverse negative blast hits (breaks bamtools if not fixed)
    if int(subjectStart) > int(subjectEnd):
        temp_subjectStart = subjectEnd
        temp_subjectEnd = subjectStart
        out_format="%s\t%s\tITS_blast_hit_%d\t%s\t%s\t.\t+\t.\tITS_blast_hits_region\n" %(subjectId,\
                            prefix, blast_count,\
                            temp_subjectStart,temp_subjectEnd)
    else:
        #direction find. ready for writing out. 
        out_format= "%s\t%s\tITS_blast_hit_%d\t%s\t%s\t.\t+\t.\tITS_blast_hits_region\n" %(subjectId,\
                                prefix, blast_count, subjectStart,\
                                subjectEnd)
    #print out_format
    return out_format


def write_out_ITS_GFF(blast, prefix, out):
    """function to write out the ITS blast hits in a GFF3
    like manner. """
    # call function to get list of blast hits.
    try:
        blast_hits = parse_blast_tab_outfile(blast)
    except:
        raise ValueError("something wrong with blast out file")
    GFF_out = open(out, "w")
    # counter to index the blast hits in the GFF file
    blast_count = 0
    # santity check to remove duplicate events
    already_seen_set = set([])

    for i in blast_hits:
        if i.startswith("#"):
            #allows the outfile to have comment lines.
            continue
        if not i.strip():
            continue #if the last line is blank
        blast_count = blast_count +1
        # check this is a unique blast hit. Remove duplicates!
        if i not in already_seen_set:
            #add this to seen set. 
            already_seen_set.add(i)
            if len(i.split("\t")) > 12:
                #remove tax id and extra coloumns - not needed.
                i = i[:12]
            if len(i.split("\t")) >12:
                raise ValueError("""custom BLAST output?
                        not enough coloumns in blast file.""")
        out_format = spit_blast_data(i, blast_count)
        #write to file
        GFF_out.write(out_format)
    #close the write file
    GFF_out.close()         


###########################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print ("v0.0.1")
    sys.exit(0)


usage = """Use as follows:
Title:
script to generate gff for ITS region BLAST hits)"
The BLAST should already have been perfomed:
 blastn -query ITS.fasta -db genome.fasta -outfmt 6 -out ITS_vs_geome.out


$ generate_ITS_GFF.py -b blast.out --prefix p.infestans -o gff.out

Note:
coloumns 6,7 and 8 are 'made' up for the purpose of this task.
BLAST hits on the negative strand will be inverted so the
start is always less than the end coordinate.

"""

parser = OptionParser(usage=usage)

parser.add_option("-b", "--blast", dest="blast", default="outfmt6.out",
                  help="the tab out file from the BLAST search",
                  metavar="FILE")
parser.add_option("--prefix", dest="prefix",
                  default="temp_name",
                  help="name for column 2 in GFF. Best to "
                  " use the origin of the data")

parser.add_option("-o", "--out_file", dest="out_file",
                  default="ITS_GFF.out",
                  help="outfile for the ITS regions in GFF format")


(options, args) = parser.parse_args()


blast = options.blast
prefix = options.prefix
out_file = options.out_file



#run the program

if not os.path.isfile(blast):
    print("sorry, couldn't open the file: " + ex.strerror + "\n")
    print ("current working directory is :", os.getcwd() + "\n")
    print ("files are :", [f for f in os.listdir('.')])
    sys_exit("\n\nInput blast file not found: %s" % blast)

# call the top function    
write_out_ITS_GFF(blast, prefix, out_file)


