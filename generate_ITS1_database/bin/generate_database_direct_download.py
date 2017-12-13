#!/usr/bin/env python
# script to return the sequences from a tab separated table of
# iTS enteries of interest.
#
# (c) The James Hutton Institute 2016-2017
# Author: Leighton Pritchard, Peter Thorpe
import os
from sys import stdin,argv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
from optparse import OptionParser
from Bio import Entrez, SeqIO
import re


def get_table_data(infile):
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
            # all the title fields in the database
            Species, date_described, marker, version, \
                     ITS_GenBank_Accession, reference, \
                     date_genbank_entry, Tax_ID, Marker_complete,\
                     genbank_entry_len, start, stop, \
                     Comments_on_ITS_region, ITS_Clade,\
                     sub_clade, UK, Geographical_distribution,\
                     Recorded_Host, Preferred_environment, \
                     DNA_only, Comment_General, \
                     Temp_Optima, Temp_Maxima, \
                     Cox, Evidence_type, Assertion_method, \
                     Source, Xref_accession, Confidence,\
                     Notes, Change_note = line.split("\t")
            line = line.rstrip()
            # rmove any spaces in name to get: Phytophthora_asiatica
            Species_name = Species.rstrip().replace(" ", "_")
            # creat dict: Phytophthora asiatica = 826
            accession_to_species[ITS_GenBank_Accession.rstrip()] = Species_name.rstrip()
            # the clade and sub clade if info is there
            clade_sub_clade = ITS_Clade + sub_clade
            species_to_clade[Species_name.rstrip()] = clade_sub_clade.rstrip()
            ITS_GenBank_Accession = ITS_GenBank_Accession.lstrip()
            accession_list.append(ITS_GenBank_Accession.rstrip())
    return accession_to_species, species_to_clade, accession_list


def harvest_desciption_for_info(descip):
    """function to split to line up at extra > sympbols.
    Returns a list split at the > symbol"""
    return descip.split(">")


def get_info(fasta_accession,
             accession_to_species,
             species_to_clade,
             genbank_species):
    """func to return info required. Takes in tow dictionaries"""
    species = accession_to_species[fasta_accession.split(".")[0]]
    # check the specie name is similar enough to the genbank name entry
    # name is collected from the GenBank download. Upper is reduce X x errors
    if species.upper() not in genbank_species.upper():
        print("%s\t%s\t%s\tDoes not match Genbank entry" % (species,
                                                            genbank_species,
                                                            fasta_accession.split(".")[0]))
    clade = species_to_clade[species]
    seq_id = "%s_%s_%s" % (clade,
                           species,
                           fasta_accession.split(".")[0])
    description = ""
    return species, clade, seq_id, description


def download_accessions(accessions):
    """func takes in a list of accessions and return
    records. """
    # Acquire the NCBI records for each accession as SeqRecord objects
    # Set your email
    Entrez.email = "peter.thorpe@hutton.ac.uk"
    idlist = ",".join(accessions)
    handle = Entrez.efetch(db="nucleotide", id=idlist, rettype="gb",
                       retmode="text")
    records = list(SeqIO.parse(handle, "genbank"))
    print("Downloaded %d GenBank records from NCBI" % len(records))
    return records


def seq_getter(tab_file, out_file):
    "function to write out the fasta file"
    # call function to populate disct
    accession_to_species, species_to_clade, accession_list\
                          = get_table_data(tab_file)
    f_out = open(out_file, 'w')
    wanted_set = set([])
    # creat a set of accessions
    records = download_accessions(accession_list)
    for entry in accession_list:
        wanted_set.add(entry.split(".")[0])
    name_set = set([])
    for seq_record in records:
        # take the first 3 element of the description. Usaully 2 are the
        # species name. But some time 3.
        genbank_species = ("_".join(seq_record.description.split()[:3]))
        fasta_accession = seq_record.id.split(".")[0]
        if fasta_accession.split(".")[0] in wanted_set:
            name_set.add(fasta_accession.split(".")[0])
            species, clade, seq_record.id, \
                     seq_record.description = get_info(fasta_accession,
                                                       accession_to_species,
                                                       species_to_clade,
                                                       genbank_species)
            seq_record.description = ""
            SeqIO.write(seq_record, f_out, "fasta")
    not_found = wanted_set.difference(name_set)
    print("#DOWNLOAD INFO")
    print("Wanted number = %d" % (len(wanted_set)))
    print("we found %d" % (len(wanted_set.intersection(name_set))))
    print("did not get %s" % not_found)
    print("\n")
    print("WANTED but not found!!:")
    for missing in not_found:
        species = accession_to_species[missing]
        clade = species_to_clade[species]
        print("%s_%s_%s" % (clade, species, missing))
    f_out.close()


usage = """Use as follows:
$ generate_database.py -t tabfile.txt -o out.fasta
requires bippython
"""

parser = OptionParser(usage=usage)

parser.add_option("-o", "--output", dest="out_file", default=None,
                  help="Output filename",
                  metavar="FILE")
parser.add_option("-t",  dest="tab_file", default=None,
                  help="tab_file containing the database informtation")


(options, args) = parser.parse_args()

out_file = options.out_file
tab_file = options.tab_file

(options, args) = parser.parse_args()

if __name__ == '__main__':
    seq_getter(tab_file, out_file)


