#!/usr/bin/env python

# title: prepare gi to desctiption database

# default 12 colounm tab output does not have description included

# author: Peter Thorpe September 2015. The James Hutton Insitute, Dundee, UK.

#i mports
import os
import sys
from optparse import OptionParser

########################################################################


def gi_to_description_generator(filename1, outfile):
    """opens up a fasta file and makes a tab separeted
    databse of gi to description for use with the diamond
    to tax info program."""
    f_out = open(outfile, "w")
    nr_fasta = open(filename1, "r")
    for line in nr_fasta:
        if line.startswith(">"):
            line = line.rstrip("\n")
            gi_number = line.split("|")[1]
            try:
                # There are some formatted description that work like this
                description = line.split("| ")[1]
            except:
                description = line.split(" ")[1:]
            try:
                data_formatted = "%s\t%s\n" %(gi_number,
                                              description.split(" >")[0])
            except:
                name = ""
                for i in description:
                    if i == ">":
                        break
                    else:
                        name = name + i
                data_formatted = "%s\t%s\n" %(gi_number, name)
            # print (data_formatted)
            f_out.write(data_formatted)
    nr_fasta.close()
    f_out.close()


#########################################################################

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.2"
    sys.exit(0)


usage = """Use as follows:

$ python prepare_gi_to_description_database.py -i nr.faa -o gi_to_des.tab

To prepare NCBI nr.faa:
blastdbcmd -entry 'all' -db nr > nr.faa

This script creates a datanse of GI number to description for the purpose
of post annotating diamond or BLAST tabular output.
Is it essential you prepare the databse before running the next script.

"""
parser = OptionParser(usage=usage)

parser.add_option("-i", "--fasta",
                  dest="nr_fasta_file",
                  default="nr.faa",
                  help="nr_fasta_file, generate using " +
                  "blastdbcmd -entry 'all' -db nr > nr.faa ,  " +
                  "you may need to: " +
                  "export BLASTDB=/PATH/TO/ncbi/extracted " +
                  "default=nr.faa")

parser.add_option("-o", "--out", dest="outfile",
                  default="gi_to_des.tab",
                  help="Output filename: default=gi_to_des.tab",
                  metavar="FILE")


# Run as script
if __name__ == '__main__':
    (options, args) = parser.parse_args()
    nr_fasta_file = options.nr_fasta_file
    outfile= options.outfile
    gi_to_description_generator(nr_fasta_file, outfile)


# more notes
"""############################################################################################
Some notes on using Diamond:


# script to get the latest NR database and NT database and make a diamond blastdatabse.


# to install diamond from source
export BLASTDB=/PATH/TO/ncbi/extracted


blastdbcmd -entry 'all' -db nr > nr.faa

/diamond-0.7.9/bin/diamond makedb --in nr.faa -d nr

diamond makedb --in uniprot_sprot.faa -d uniprot

diamond makedb --in uniref90.faa -d uniref90

covert output to tab:
$ diamond view -a diamond.daa -f tab -o name.tab


from stdin:
diamond makedb --in /dev/stdin -d tiny_from_stdin < tiny.faa
"""
