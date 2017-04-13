#!/usr/bin/env python
# title: prepare gi to desctiption database
# default 12 colounm tab output does not have description included
# author: Peter Thorpe September 2015. The James Hutton Insitute, Dundee, UK.

# imports
import os
import sys
from optparse import OptionParser

########################################################################


def get_annotation(line, descriptions=4):
    """function to get the associated annotation with its GI number
    takes in a line and return the number of 'descriptions'
    annotations.
    The middle of the code is messy to deal with varying annotation
    methods found in the nr.faa file.
    returns a string of space sep annotations.
    """
    number_of_descriptions = int(descriptions)
    line = line.split(">")
    annot = ""
    count = 0
    # get the number of descriptions in the fa header
    # remeber elemnt 1 is blank
    number_of_hits = len(line) - 1
    if number_of_hits < descriptions:
        # if this is less than the amount we want, have to
        # adjust or will get None back
        descriptions = number_of_hits
    for entry in line:
        if entry == "":
            continue
        # print ("entry = ", entry)
        if "| " in entry:
            # e.g.>gi|489223532|ref|WP_003131952.1| 30S ribosomal protein
            description = entry.split("| ")[1]
        elif ": " in entry:
            # e.g.>gi|489223532|ref|WP_003131952.1: annotation
            description = entry.split(": ")[1]
        else:
            # SOME ARE FORMATTED JUST WITH A SPACE
            # e.g.>gi|489223532|ref|WP_003131952.1 annotation
            description = " ".join(entry.split(" ")[1:])
        if count == 0:
            annot = description  # first time so we dont start with a space
        else:
            # append the new annot
            annot = annot + " " + description.rstrip()
        # this will get us the first 4 blast descriptions. default
        # can be modified
        count = count + 1
        if count == descriptions:
            return annot  # string


def get_accession_number(line):
    """function to get the GI number in the fa line.
    take in the line, returns the GI number of the first fa entry
    Example:
    >gi|446057344|ref|WP_000135199.1| MULTISPECIES:
    30S ribosomal protein S18"""
    acc_number = line.split("| ")[0]
    return acc_number.split("|")[3]


def acc_to_description_generator(filename1, descriptions,
                                 outfile):
    """opens up a fasta file and makes a tab separeted
    databse of gi to description for use with the diamond
    to tax info program."""
    f_out = open(outfile, "w")
    nr_fasta = open(filename1, "r")
    for line in nr_fasta:
        if line.startswith("#"):
            continue  # comment line
        if not line.strip():
            continue  # if the last line is blank
        if line.startswith(">"):
            line = line.rstrip("\n")
            acc_number = get_accession_number(line)
            annot = get_annotation(line, descriptions)
            try:
                data_formatted = "%s\t%s\n" %(acc_number,
                                              annot)
            except ValueError:
                print ("something failed in getting the descriptions")
                os.exit()
            # print (data_formatted)
            f_out.write(data_formatted)
    nr_fasta.close()
    f_out.close()


#########################################################################

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.3"
    sys.exit(0)


usage = """Use as follows:

$ python prepare_accession_to_description_db.py -i nr.faa -o acc_to_des.tab

To prepare NCBI nr.faa:
blastdbcmd -entry 'all' -db nr > nr.faa

This script creates a database of accession number to description for the purpose
of post annotating diamond or BLAST tabular output.
Is it essential you prepare the databse before running the next script.

Use the --descriptions options to control how many descriptions for the same hit
are returned. default 4 hits.


descriptions = number of blast descriptions (default is 4) to return which are asscosiated with this accession number. 
More will be useful but your
file will get very big. Quickly!
./prepare_accession_to_description_db.py -h
Usage: Use as follows:

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

parser.add_option("-n", "--descriptions", dest="descriptions",
                  default=4,
                  help="number of blast descriptions to return " +
                  " which are asscosiated with this GI number")

parser.add_option("-o", "--out", dest="outfile",
                  default="acc_to_des.tab",
                  help="Output filename: default=acc_to_des.tab",
                  metavar="FILE")


# Run as script
if __name__ == '__main__':
    (options, args) = parser.parse_args()
    nr_fasta_file = options.nr_fasta_file
    if not os.path.isfile(nr_fasta_file):
        print ("sorry cannot find you %s file" % nr_fasta_file)
        print ("please check this command again, " +
               "with the full path if required")
        os._exit(0)
    outfile = options.outfile
    descriptions = options.descriptions
    acc_to_description_generator(nr_fasta_file, descriptions, outfile)


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
