#!/usr/bin/env python

# title: script to pull out the regions/ domains of interest from a
# hmmsearch --domtblout - NUCLEOTIDE VERSION!!!!!!!!!
# NOTE THIS IS AN OLD SCRIPT THAT NEEDS UPDATING

# imports
# biopython
# os imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#################################################################

def parse_file(HMM_search_file):
    """function to open and parse a file.
    returns """
    f = open(HMM_search_file, "r")
    # assign the file contents to the variable data
    data = f.readlines()
    # remove the \n new line and \t characters
    data1 = [line.rstrip("\n").split() for line in (data)
             if line.strip() != "" and not line.startswith("#")]
    HMM_search_data = [(str(s[0]), int(s[17]), \
                        int(s[18]), int(s[2])) \
                       for s in (data1)]
    f.close()
    return HMM_search_data


def domain_getter(filename, HMM_search_file, outfile):
    """function to write the sequence identifed in the
    domain regions of a pfam domain to a file
    - required hmmsearch --domtblout -
    NUCLEOTIDE VERSION"""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    domain_tbl_out = open("seqid_domain_coordinates.txt", "w")
    # open to write to stat stop coordinate to
    HMM_search_data = parse_file(HMM_search_file)
    f_out= open(outfile, 'w')
    for seq_record in SeqIO.parse(filename, "fasta"):
        for i in HMM_search_data:
            HMM_search_name = i[0]
            # ali coord from and to from the HMM out file is what we want
            HMM_search_position_start = ((i[1]))
            HMM_search_position_stop = (i[2])
            HMM_search_position_start_real = HMM_search_position_start
            seq_length = i[3]
            #print HMM_search_name
            if HMM_search_name == seq_record.id:
                assert HMM_search_position_start_real \
                       < HMM_search_position_stop <= len(seq_record), \
                       "HMM_name %s, Record %s lenh %i, coords %i to %i" \
                       % (HMM_search_name, seq_record.id,
                          len(seq_record),\
                          HMM_search_position_start_real,
                          HMM_search_position_stop)
                #if seq_length == len(seq_record):
                #print seq_record.id
                temp_name = " ".join(seq_record.id.split("_")[1:-1])
                domain = "%s\t%s\t%s\n" % (temp_name,
                                           HMM_search_position_start,
                                           HMM_search_position_stop)
                domain_tbl_out.write(domain)
                ITS = seq_record.seq[HMM_search_position_start_real:HMM_search_position_stop]

                output_formatted = '>%s\n%s\n' %(seq_record.id,
                                                 ITS)
                f_out.write(output_formatted)
    f_out.close()
    domain_tbl_out.close()


##########################################################################

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)

usage = """Use as follows:

Script to get the NUCLEOTIDE domains/ region of interest as
defined in the hmmsearch --domtblout

$ python get_DOMAIN_region_I_want_from_fasta_Nucleotide.py
-i in.fasta --hmm ITS1_hmm_domain_table.out -o outfile_domains.fasta

requires bippython
"""

parser = OptionParser(usage=usage)

parser.add_option("-i", "--in", dest="in_file",
                  default=None,
                  help="nt_file used to generate the AA file"
                  "used for the hmmsearch")
parser.add_option("-o", "--output", dest="out_file",
                  default="nt_domains.fasta",
                  help="Output filename - domains only",
                  metavar="FILE")
parser.add_option("--hmm", dest="hmm_output_file",
                  default=None,
                  help="hmm_output_file --domtblout  output")


(options, args) = parser.parse_args()

filename = options.in_file
hmm_output_file = options.hmm_output_file
out_file = options.out_file

if __name__ == '__main__':
    domain_getter(filename, hmm_output_file, out_file)






