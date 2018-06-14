#!/usr/bin/env python
# title: Populate the Metapy results for biologists to look at it in excel
# (c) The James Hutton Institute 2018
# Author: Peter Thorpe, Leighton Pritchard
import sys
import os
import argparse


VERSION = "summerise results: v0.01"
if "--version" in sys.argv:
    print(VERSION)
    sys.exit(1)

usage = """

Save the lists of sample with the hosts into a text file:
#plate_well	Sample_Number	Sample_type_(root/plant/water/blank/other)	Host_Common_Name
2D	N010-161128-2A-F_S38_L001	Filter		Reservoir used for whole site before NaOCl treatment

This is an example. This will go through all the results in a folder
for a tool of your choice and populate a new text file with
the old info and add the species found. But only those with more reads
that the given threshold. Dafault 50.

 python populate_excel_sheet_swarm.py -t 100 -i infile.tx -o outfile.text
 
"""
if "--help" or "-h" in sys.argv:
    print(usage)


def get_args():
    parser = argparse.ArgumentParser(description="Pipeline: cluster " +
                                     "data for metabarcoding ",
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("metapy")[0]
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-t", "--threshold", dest='threshold',
                          action="store", default=10,
                          type=int,
                          help="min number of reads to consider this an" +
                          " actual hit. Default 50")

    optional.add_argument("-i", "--in", dest='infile',
                          action="store",
                          default="Illumina01112017.txt",
                          type=str,
                          help="the tab file to fill in. Input file")

    optional.add_argument("-p", "--program", dest='program',
                          action="store",
                          default="bowtie",
                          type=str,
                          help="program of interest to get results for" +
                          ". Default bowtie.")

    optional.add_argument("-o", "--out", dest='out',
                          action="store",
                          default="Illumina_results.txt",
                          type=str,
                          help="the tab file to fill in. Input file")

    optional.add_argument('--version',
                          action='version',
                          version="%s: metapy.py " + VERSION)
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
    # used to be [1] for thapbi project
    sample_name = data[0]
    sample_set.add(sample_name.rstrip())  # all kinds of whitespace
    return sample_set


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


def parse_swarm_result(swarm_result_file, threshold):
    """funct to parse the swarm result file and return
    a tab separated list of species found .
    infile:     # cluster_number species number_of_reads_hitting_species
    returns a string of species\tnum_reads etc ...
    threshold = min number of reads clustering with the phy to be considered
    an actual hit
    """
    result = ""
    with open(swarm_result_file) as handle:
        for line in handle:
            if test_line(line):
                cluster_number, species, number_of_reads = line.rstrip().split("\t")
                cluster_result = "%s\t%s\t" % (species.strip(), number_of_reads.strip())
                number_of_reads = number_of_reads.strip()
                if int(number_of_reads) >= threshold:
                    result = result + cluster_result
    return result


def write_out_result(indata, outfile):
    """takes in string, writres out to file"""
    f_out = open(outfile, "w")
    title = "#plate_well\tSample_Number\tSample_type_(root/plant/water/blank/other)" +\
            "\tHost_Common_Name\tHost_latin_name/Sampling_Name\tPhytophthora_species\t" +\
            "num_reads\tPhytophthora_species\tnum_reads\n"
    #f_out.write(title)
    #indata = indata.replace("_abundance=1", "")
    print "indata = ", indata
    for result in indata:
        f_out.write(result)
    f_out.close()


def populate_result_list(full_data, entry, result,
                         full_data_with_phytophora):
    """takes in the original file already converted to a list.
    Take in the phy result string.
    resturns list to write to file in other function
    """
    for data in full_data:
        data = data.rstrip("\n")
        # used to be [1] fr thapbi project
        if len(data.split("\t"))> 1:
            data_list = data.split("\t")
            sample = data_list[0]
        else:
            sample = data
        if sample.strip().replace("_RESULTS", "") == entry.strip().replace("_bowtie_perfect.RESULTS", ""):
            #print "SAMPLE:", sample, "Entry", entry, "results", result

            full_data_with_phytophora = sample + "\t" + full_data_with_phytophora + \
                                        data.strip() + "\t" + result.rstrip() + "\n"
            if "PCR" in sample:
                print "full_data_with_phytophora = ", full_data_with_phytophora
    return full_data_with_phytophora


args, FILE_DIRECTORY = get_args()
# Run as script
if __name__ == '__main__':
    PROG_OF_INTEREST = args.program
    # call the function to get a list of results wanted
    full_data, sample_set = parse_text_file(args.infile)
    directory = "."
    full_data_with_phytophora = ""
    # get all folder names os.walk(directory)
    results_folders = [x[0] for x in os.walk(directory)]
    for result in results_folders:
        if PROG_OF_INTEREST in result:
            if "DADA2" not in result:
                for entry in sample_set:
                    if entry.endswith("_RESULTS"):
                        entry = entry.split("_RESULTS")[0]
                    wanted_file = entry + "_bowtie_perfect.RESULTS"
                    #print wanted_file
                    wanted_file_full = os.path.join(result, wanted_file)
                    if os.path.isfile(wanted_file_full):
                        # print wanted_file_full
                        # result is string of species\tnum_reads etc ...
                        result = parse_swarm_result(wanted_file_full, args.threshold)
                        # full data is the original file converted to a list
                        # entry is the sample of interets.
                        # result in the phy spcies and num of reads.
                        full_data_with_phytophora = populate_result_list(full_data,
                                                                         entry,
                                                                         result,
                                                                         full_data_with_phytophora)
    # call the function to write out the full_data_with_phytophora
    write_out_result(full_data_with_phytophora, args.out)


