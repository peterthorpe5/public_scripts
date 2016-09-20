#!/usr/bin/env python

######################################################################################################################
#Title: pick longest transdecoder component
######################################################################################################################
"""
why? Trandecoder can predict multiple cds per transcript.
This script pick the longest

"""

#######################################################################################################################
# imports for system

import sys
import os
from optparse import OptionParser

##########################################################################################################################
#imports for functions
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import tempfile
from collections import deque
import datetime

#make a temp_folder_for_all_the_out_files
file_name = 'test.txt'
working_dir = os.getcwd()
dest_dir = os.path.join(working_dir, 'longest_representative')
try:
    os.makedirs(dest_dir)
except OSError:
    print ("folder already exists, I will write over what is in there!!")

#########################################################################################################################
# functions 
##########################################################################################################################

def sys_exit(msg, error_level=1):
   """Print error message to stdout and quit with given error level."""
   sys.stderr.write("%s\n" % msg)
   sys.exit(error_level)


def strip_to_match_transcript_name(identifier):
    "remove the cds and split at the first pipe"
    return identifier.replace("cds.", "").split("|m.",1)[0]

def find_longest_components(filename1, cds_database, out_filename):
    """this is a function to open up a fasta file, and
    producea a list of the longest representative transcripts per gene"""

    #this is out list of so called longest matches which we will append and remove as applicable
    top_hits = []
    #current sequence lenth score value "to beat"
    current_lenth = int(0)
    #set up variables to assgn lastest values to ...
    transcriptome_Genes_names = set([])
    last_gene_name = ""
    last_component = ""
    loop_count = 0
    for seq_record in SeqIO.parse(filename1, "fasta"):
        sequence_len = len(seq_record)
        sequence_name = seq_record.id
        component = strip_to_match_transcript_name(sequence_name)
        #first time we see any record, save the values:
        if loop_count == 0:
            loop_count = loop_count+1
            last_gene_name = sequence_name
            current_lenth = sequence_len
            last_component= component
            top_hits.append(seq_record.id)
            
        ##############################################################################
        #first block: if the names are the same, is the new length of sequence longer?
        if component == last_component:
            #print "yes:", component, "component",  last_component, "last_component", seq_record.id
            #print "current_lenth", current_lenth
            if sequence_len > current_lenth:
                #print "sequence_len > current_lenth", sequence_len, current_lenth
                del top_hits[-1]
                top_hits.append(seq_record.id) 
        #############################################################################
        # second block: if the name is new, put it in the name set.
        # use this sequence-length as the new one to "beat"
        else:
            top_hits.append(seq_record.id)
            last_gene_name = sequence_name
            current_lenth = sequence_len
            last_component= component
    outfile = open(out_filename,"w")
    for i in top_hits:
        seq_record =  cds_database[i]
        SeqIO.write(seq_record, outfile, "fasta")
    outfile.close()
    cds_database_new = SeqIO.index(out_filename, "fasta", key_function=strip_to_match_transcript_name)
    return cds_database_new

def parse_predicted_CDS_file(cds_file, out):
    """parse the cds file and index it"""
    # this is for transdecoder names. May need to alter for other tools
    try:
        cds_database = SeqIO.index(cds_file, "fasta", key_function=strip_to_match_transcript_name)
        return cds_database
    except ValueError:
        print ("""looks like multiple cds were predicted per transcript - cannot change names.
        I am going to pick the longest representative cds per transcripts. I only do this if there are
        multiple cds predicted per transcript, otherwise this message is not shown""")
        cds_database = SeqIO.index(cds_file, "fasta")
        #basically there are duplicates for each transcript. So, find the longest
        #representative and create a new cds_database, based on that
    #call function
    cds_database_new = find_longest_components(cds_file, cds_database, "./longest_representative/%s" % out)
    #return a seq_record object that can be accessed in a dictionary like manner
    return cds_database_new

           

###############################################################################################
#to run it:
usage = """Use as follows:

This script chooses the longest cds component per transcript
predicted by transdecoder.

python pick_longest_cds.py --cds nt_coding_seq
 results will be in ./longest_representative/out.fasta

require: Biopython

"""
parser = OptionParser(usage=usage)


parser.add_option("--cds", dest="cds", default=None,
                  help="the predicted cds from the transcriptome assembly .fasta",
                  metavar="FILE")
parser.add_option("-o", "--out", dest="outfile", default="results.out",
                  help="Output filename (default: results.out)",
                  metavar="FILE")

(options, args) = parser.parse_args()


#-o
outfile= options.outfile
#--cds
cds_file = options.cds

parse_predicted_CDS_file(cds_file, outfile)


