#!/usr/bin/env python
# title: Populate the Metapy results for biologists to look at it in excel
# (c) The James Hutton Institute 2018
# Author: Peter Thorpe, Leighton Pritchard
import sys
import os
import argparse
from collections import defaultdict


VERSION = "summerise results: v0.01"
if "--version" in sys.argv:
    print(VERSION)
    sys.exit(1)

usage = """

Save the lists of sample with the hosts into a text file:
#plate_well	Sample_Number	Sample_type_(root/plant/water/blank/other)	Host_Common_Name
2D	N010-161128-2A-F_S38_L001	Filter		Reservoir used for whole site before NaOCl treatment

This is an example. This will go through all the "novel blast " results
in a folder  and populate a new text file with
the old info and add the species found. But only those with more reads
that the given threshold. Default 5 mismatches (-m)

 python populate_excel_sheet_with_novel_blast_results.py -m 5 -i infile.tx -o outfile.text

"""
#if "--help" or "-h" in sys.argv:
    #print(usage)


def get_args():
    parser = argparse.ArgumentParser(description="Pipeline: cluster " +
                                     "data for metabarcoding ",
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("metapy")[0]
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument("-i", "--in", dest='infile',
                          action="store",
                          default="Illumina01112017.txt",
                          type=str,
                          help="the tab file to fill in. Input file")

    optional.add_argument("-p", "--program", dest='program',
                          action="store",
                          default="Swarm",
                          type=str,
                          help="program of interest to get results for" +
                          ". Default Swarm.")

    optional.add_argument("-m", "--mismatches", dest='mismatches',
                          action="store",
                          default=5,
                          type=int,
                          help="number of mismatches to filter results with")

    optional.add_argument("-o", "--out", dest='out',
                          action="store",
                          default="Illumina01112017_novel_blast_results.txt",
                          type=str,
                          help="the tab file to fill in. Input file")

    optional.add_argument('--version',
                          action='version',
                          version="%s: populate.py " + VERSION)
    args = parser.parse_args()
    return args, file_directory


def test_line(line):
    """returns true lines. Not comments or blank line"""
    if not line.strip():
        return False  # if the last line is blank
    if line.startswith("#"):
        return False  # comment line
    if line.startswith("    # "):  # swarm result file
        return False  # comment line
    if line.startswith("		p"):
        return False  # comment line
    return line.rstrip()


def split_line_return_sample(line, sample_set):
    """func to split up a line and return a list of samples
    - because result havd come with a mixture of underscores and
    hyphens in all possible plces, I will just replace these."""
    line = line.rstrip()
    data = line.split("\t")  # inconsitent data entry!! cannot assign
    # was [1] for thapbi data
    sample_name = data[0]
    sample_set.add(sample_name.rstrip())  # all kinds of whitespace
    return sample_set


def get_blast_data_element(line, novel_hit, sample_name,
                           sample_name_to_blast_hit,
                           cluster_size,
                           mismatches_threshold):
    """function to eturn the individual element of the blast output.
    #Query	mismatches	evalue	percent_identity	database
    Those which are under mismatches"""
    Query, mismatches, evalue, percent_identity, database = line.rstrip().split("\t")
    if int(mismatches) <= mismatches_threshold:
        accession_number = database.split("|")[3]
        database  = database.split("|")[4]
        species_name = " ".join(database.split()[:4])
        data = "%s%s|%s\t%s\t" % (novel_hit, str(accession_number), str(species_name), cluster_size)
        sample_name_to_blast_hit[sample_name] += data
    return sample_name_to_blast_hit


def split_line_blast_file(filename, mismatches,
                          sample_name_to_blast_hit,
                          sample_name_to_cluster_size,
                          PROGRAM="Swarm"):
    """func to split up a line and return a list of samples
    - because result havd come with a mixture of underscores and
    hyphens in all possible plces, I will just replace these."""
    sample_name = os.path.split(filename)[-1]
    cluster_size = sample_name.split("len")[1]
    cluster_size = cluster_size.split("_vs")[0]
    sample_name = sample_name.split("_" + PROGRAM)[0]
    novel_hit = ""
    cluster_string = cluster_size + ", "
    sample_name_to_cluster_size[sample_name] += cluster_string
    with open(filename) as handle:
        for line in handle:
            if test_line(line):
                line = line.rstrip()
                sample_name_to_blast_hit = get_blast_data_element(line, novel_hit,
                                                                  sample_name,
                                                                  sample_name_to_blast_hit,
                                                                  cluster_size,
                                                                  mismatches)
    return sample_name_to_blast_hit, sample_name_to_cluster_size



def parse_text_file(text_file):
    """func to open and parse the text file"""
    sample_set = set([])
    full_data = []
    with open(text_file) as handle:
        for line in handle:
            full_data.append(line)
            if test_line(line):
                sample_set = split_line_return_sample(line, sample_set)
    return full_data, sample_set


def write_out_result(indata, outfile):
    """takes in string, writres out to file"""
    f_out = open(outfile, "w")
    title = "#plate_well\tSample_Number\tSample_type_(root/plant/water/blank/other" +\
            "\tHost_Common_Name\tHost_latin_name/Sampling_Name\tcluster_sizes" +\
            "\tspecies\tcluster_size\tspecies\tcluster_size\tspecies\tcluster_size ...\n"
    f_out.write(title)
    indata = indata.replace("_abundance=1", "")
    for result in indata:
        f_out.write(result)
    f_out.close()


def populate_result_list(full_data,
                         sample_name_to_blast_hit,
                         sample_name_to_cluster_size,
                         full_data_with_phytophora):
    """takes in the original file already converted to a list.
    Take in the phy result string.
    resturns list to write to file in other function
    """
    for line in full_data:
        data = line.rstrip("\n")
        data_list = data.split("\t")
        while len(data_list) < 6:
            data_list.append("\t")
            # was [1] for thpabi data
        sample = data_list[0]
        sample = sample.replace("_RESULTS", "")
        print sample
        try:
            blast_result = sample_name_to_blast_hit[sample.rstrip()]
        except:
            blast_result = ""
        cluster_size = sample_name_to_cluster_size[sample.rstrip()]
        full_data_with_phytophora = full_data_with_phytophora + \
                                    line.rstrip() + "\t" +\
                                    cluster_size.rstrip(", ") + "\t" + \
                                    blast_result + "\n"
    return full_data_with_phytophora


args, FILE_DIRECTORY = get_args()
# Run as script
if __name__ == '__main__':
    sample_name_to_blast_hit = defaultdict(str)
    sample_name_to_cluster_size = defaultdict(str)
    PROG_OF_INTEREST = args.program
    # call the function to get a list of results wanted
    full_data, sample_set = parse_text_file(args.infile)
    directory = "."
    full_data_with_phytophora = ""
    last_name = ""
    last_hit = ""
    result = ""
    for filename in os.listdir(".") :
        if not filename.endswith(".result.txt"):
            continue
        sample_name_to_blast_hit, \
            sample_name_to_cluster_size = split_line_blast_file(filename,
                                                                args.mismatches,
                                                                sample_name_to_blast_hit,
                                                                sample_name_to_cluster_size,
                                                                PROG_OF_INTEREST)

    full_data_with_phytophora = populate_result_list(full_data,
                                                     sample_name_to_blast_hit,
                                                     sample_name_to_cluster_size,
                                                     full_data_with_phytophora)

        
        
    # call the function to write out the full_data_with_phytophora
    write_out_result(full_data_with_phytophora, args.out)


