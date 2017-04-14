#!/usr/bin/env python
# Title: Get nt_seq of interest, correspoding protein
# seq effector blast matches
# why: script to compare generated gene models against a set of
# known sequence
# author: Peter Thorpe September 2015. The James Hutton Insitute, Dundee, UK.

# imports
from __future__ import print_function
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import sys
from optparse import OptionParser


def make_folder(folder_name="known_fa_all_hits"):
    """function to make a folder with a given name.
    Take in the name, returns nothing"""
    dest_dir = os.path.join(os.getcwd(), 'known_fa_all_hits')
    try:
        os.makedirs(dest_dir)
    except OSError:
        print("already exists")


def test_line(line):
    """returns true lines. Not comments or blank line"""
    if not line.strip():
        return False  # if the last line is blank
    if line.startswith("#"):
        return False  # comment line
    return line


def open_blast_file(blast_file):
    """funtion to open the blast file.
    returns a \n separeted list
    read in the tab file. Reads whole file into memory.
    Could be altered for more efficiency
    """
    with open(blast_file) as file:
        return file.read().split("\n")


def parse_blast_tab_hit(blast_line):
    """takes in a tab blast line. Splits at \t
    retunrs the names of the query and the db
    hit"""
    data = blast_line.rstrip("\n").split("\t")
    gene = data[0]
    blast_hit_matches = data[1]
    return gene, blast_hit_matches


def get_known_name(KNOWN_FA):
    """takes in a fasta of known seq that we want to compare the
    gene models against. This funtion returns a list of names
    and a seq_IO_index of this file"""
    known_name_list = []
    print("Indexing...%s" % KNOWN_FA)
    known_seq_db = SeqIO.index(KNOWN_FA, "fasta")
    for seq_record in known_seq_db:
        known_name_list.append(seq_record)
    return known_seq_db, known_name_list


def generate_dict_of_files(known_name_list,
                           folder_name,
                           prefix):
    """funtion to create a dictionary of outfie.
    Take in a list, foldername and prefix.
    returns:
    folder_name, gene + "_PREFIX.fasta.
    All these files are also open for writing to"""
    phandles = dict()
    for gene in known_name_list:
        suffix = "_%s.fasta" % prefix
        filename = os.path.join(folder_name, gene + suffix)
        phandles[gene] = open(filename, "w")
    return phandles


# main function
def seq_getter(blast_file,
               cds_file,
               protein_file,
               known_seq_db,
               known_name_list,
               folder_name):
    """Function to convert a top blast anlaysis query versus seq
    to a file containing these sequences."""
    print("current cds_file is false")
    cds_file = False
    if cds_file:
        # if we have a nt file... the creat files associated
        # with what is in the file. + get a filename dict
        nhandles = generate_dict_of_files(known_name_list,
                                          folder_name,
                                          "cds")
    if protein_file:
        # if we have a AA file... the creat files associated
        # with what is in the file. + get a filename dict
        phandles = generate_dict_of_files(known_name_list,
                                          folder_name,
                                          "pep")
    if protein_file:
        # index
        protein_sequences = SeqIO.index(protein_file, "fasta")
    if cds_file:
        nucleotide_sequences =  SeqIO.index(cds_file, "fasta")
    print("Starting output...")
    names_already_printed = set([])
    blast_hits_wanted = open_blast_file(blast_file)
    for line in blast_hits_wanted:
        if not test_line(line):
            continue
        gene, blast_hit_matches = parse_blast_tab_hit(line)
        if cds_file:
            # get the nt seq from the gene models file
            seq_record = nucleotide_sequences[blast_hit_matches]
            SeqIO.write(seq_record, nhandles[gene], "fasta")
        if protein_file:
            # get the AA seq from the gene models file
            seq_record = protein_sequences[blast_hit_matches]
            SeqIO.write(seq_record, phandles[gene], "fasta")
        if not gene in names_already_printed:
            # get the seq from the known db
            seq_record = known_seq_db[gene]
            SeqIO.write(seq_record, phandles[gene], "fasta")
            names_already_printed.add(gene)
    # loop to close all the open files. There could be many!!
    if cds_file:
        for gene in known_name_list:
            nhandles[gene].close()
    for gene in known_name_list:
        phandles[gene].close()

#################################################################################################

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.2")
    sys.exit(0)


usage = """Use as follows:

python ${python_directory}/Get_sequence_from_tab_blast.py
	    -b test_fa_vs_known_fa.tab
	    --known_fa database.fa
	    -p aa.fasta
	    -n nt.fa

Script to creat a separted file for all the top blast matches.
Why? So they can be aligned an get a visual representation of
how good the gene models are
"""

parser = OptionParser(usage=usage)

parser.add_option("-b", dest="blast_output",
                  default=None,
                  help="in xml")

parser.add_option("--known_fa", dest="known_fa",
                  default=None,
                  help="this is the seq database of known seq " +
                  "to compare to ",
                  metavar="FILE")

parser.add_option("-p", dest="protein_file",
                  default=None,
                  help="predicted amino acid seq of genes")

parser.add_option("-f","--folder",
                  dest="folder_out",
                  default="known_fa_all_hits",
                  help="folder where sequences will be dumped")

parser.add_option("-n", dest="nuc",
                  default=None,
                  help="predcited nt cds of genes. ",
                  metavar="FILE")


# Run as script
if __name__ == '__main__':
    (options, args) = parser.parse_args()
    KNOWN_FA = options.known_fa
    BLAST_OUTPUT = options.blast_output
    FOLDER_OUT = options.folder_out
    PROTEIN_FILE = options.protein_file
    NT_FILE = options.nuc

    # call functions
    make_folder(FOLDER_OUT)
    known_seq_db, known_name_list = get_known_name(KNOWN_FA)
    seq_getter(BLAST_OUTPUT,
               NT_FILE,
               PROTEIN_FILE,
               known_seq_db,
               known_name_list,
               FOLDER_OUT)
