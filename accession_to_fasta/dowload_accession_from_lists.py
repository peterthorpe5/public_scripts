#!/usr/bin/env python
# script to return the sequences from a list of accession
#
# (c) The James Hutton Institute 2016-2017
# Author: Peter Thorpe
import os
from sys import stdin,argv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
from optparse import OptionParser
from Bio import Entrez, SeqIO


def get_accession_list(infile):
    """function to open a tab sep table and return
    info required.  Returns dictionaries.
    species_to_clade  accession_to_species"""
    accession_to_species = dict()
    species_to_clade = dict()
    accession_list = []
    with open(infile) as handle:
        for line in handle:
            # print(line)
            if line.startswith("#"):
                continue
            if not line.strip():
                continue # if the last line is blank
            line = line.rstrip()
            accession_list.append(line)
    print("loaded accssions: eg.", accession_list[:2], " ...")
    return accession_list



def download_accessions(accessions, email):
    """func takes in a list of accessions and return
    records. """
    # Acquire the NCBI records for each accession as SeqRecord objects
    # Set your email
    print("Dowloading data ...this may take a while")
    
    Entrez.email = str(email)
    idlist = ",".join(accessions)
    handle = Entrez.efetch(db="nucleotide", id=idlist, rettype="gb",
                       retmode="text")
    records = list(SeqIO.parse(handle, "genbank"))
    print("Downloaded %d GenBank records from NCBI" % len(records))
    return records


def seq_getter(tab_file, email, out_file):
    "function to write out the fasta file"
    # call function to populate disct
    accession_list = get_accession_list(tab_file)
    f_out = open(out_file, 'w')
    wanted_set = set([])
    # creat a set of accessions
    records = download_accessions(accession_list, email)
    for entry in accession_list:
        wanted_set.add(entry)
    name_set = set([])
    for seq_record in records:
        # take the first 3 element of the description. Usaully 2 are the
        # species name. But some time 3.
        genbank_species = ("_".join(seq_record.description.split()[:3]))
        fasta_accession = seq_record.id.split(".")[0]
        # print("genbank_species = \t%s\t%s" % (genbank_species, fasta_accession))
        name_set.add(fasta_accession.split(".")[0])
        seq_record.description = ""
        SeqIO.write(seq_record, f_out, "fasta")
    not_found = wanted_set.difference(name_set)
    print("#DOWNLOAD INFO")
    print("Wanted number = %d" % (len(wanted_set)))
    print("we found %d" % (len(wanted_set.intersection(name_set))))
    print("did not get %s" % len(not_found))
    print("\n")
    print("WANTED but not found!!:")
    f_out.close()


usage = """Use as follows:
$ generate_database.py -t tabfile.txt -o out.fasta
requires bippython

tabfile, is simply a list of accession.

one per line
"""

parser = OptionParser(usage=usage)

parser.add_option("-o", "--output", dest="out_file",
                  default="sequences.fasta,
                  help="Output filename",
                  metavar="FILE")
parser.add_option("-e",  dest="email",
                  default=None,
                  help="your email address, this is essential for NCBI")
parser.add_option("-t",  dest="tab_file",
                  default=None,
                  help="list of accession of interest")


(options, args) = parser.parse_args()

out_file = options.out_file
tab_file = options.tab_file
email = options.email

(options, args) = parser.parse_args()

if __name__ == '__main__':
    seq_getter(tab_file, email, out_file)


