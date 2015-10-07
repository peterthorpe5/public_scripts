#!/usr/bin/env python

#title: script to pull out the regions/ domains of interest from a
#hmmsearch --domtblout - NUCLEOTIDE VERSION!!!!!!!!!

#imports
#biopython
#os imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser

###############################################################################
def domain_getter(filename, HMM_search_file, outfile):
    """"this function print the domain regions of a pfam domain to a file
    - required hmmsearch --domtblout - NUCLEOTIDE VERSION"""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    f= open(HMM_search_file, "r")
    #assign the file contents to the variable data
    data = f.readlines()
    #remove the \n new line and \t characters
    data1 = [line.rstrip("\n").split() for line in (data)
             if line.strip() != "" and not line.startswith("#")]
    #to convert data 1 into a list of tuples.
    #remove the title of the file
    #data2 = data1 [:]
    #print data2
    #THE NEXT LINE IS SPECIFIC TO THE OVERAL TASK NOT TO THIS FUNCTION
    HMM_search_data = [(str(s[0]), int(s[17]), int(s[18]),int(s[2])) for s in (data1)]
    #print HMM_search_data
    from Bio import SeqIO
    f= open(outfile, 'w')

        #print seq_record.id
        #THIS get the names from HMM_search file
    for seq_record in SeqIO.parse(filename, "fasta"):
        for i in HMM_search_data:
            HMM_search_name = i[0]
            HMM_search_position_start = (3*(i[1]))-3
            HMM_search_position_stop = 3*(i[2])
            HMM_search_position_start_real = HMM_search_position_start
            seq_length = i[3]
            #print HMM_search_name
            if HMM_search_name == seq_record.id:
                assert HMM_search_position_start_real < HMM_search_position_stop <= len(seq_record), \
                       "HMM_searchname %s, Record %s length %i, coords %i to %i" \
                       % (HMM_search_name, seq_record.id, len(seq_record),\
                          HMM_search_position_start_real, HMM_search_position_stop)
                if seq_length == len(seq_record):

                    output_formatted = '>%s\t%i:%i\n%s\n' %(seq_record.id, HMM_search_position_start,\
                                                       HMM_search_position_stop,\
                                    seq_record.seq[HMM_search_position_start_real:HMM_search_position_stop])

                    f.write(output_formatted)

    f.close()
    return True


#################################################################################################################################################

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

Script to get the NUCLEOTIDE domains/ region of interest as defined in the hmmsearch --domtblout

$ python get_DOMAIN_region_I_want_from_fasta_amino_acid.py -i in.fasta -o outfile_domains.fatsa

requires bippython
"""

parser = OptionParser(usage=usage)

parser.add_option("-i", "--in", dest="in_file", default=None,
                  help="nt_file used to generate the AA file used for the hmmsearch")
parser.add_option("-o", "--output", dest="out_file", default=None,
                  help="Output filename - domains only",
                  metavar="FILE")
parser.add_option("--hmm", dest="hmm_output_file", default=None,
                  help="hmm_output_file --domtblout  output")




(options, args) = parser.parse_args()

filename = options.in_file
hmm_output_file = options.hmm_output_file
out_file = options.out_file

(options, args) = parser.parse_args()



domain_getter(filename, hmm_output_file, out_file)

print 'done'





