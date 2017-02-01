#title: Parse clusters and find the phy spcies
#in the cluster
#author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

#imports
import os
import sys
from optparse import OptionParser
import datetime
import os
from sys import stdin,argv

import gzip
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator



#####################################################################

def turn_fq_to_dic(in_fastq, barcode_length):
    """function to turn fq file into a dict
    read = seq[:length of bar code]
    return dict """
    #THIS FUNCTION IS NOW NOT USED.
    read_to_barcode_dict = dict()
    #open the fastq file
    in_file = open(in_fastq)
    # iterate through the fatsq file
    total_reads = 0
    for (title, sequence, quality) in FastqGeneralIterator(in_file):
        total_reads = total_reads+1
        title = title.split(" ")[0]
        read_to_barcode_dict[title.split(" ")[0]] = sequence[:barcode_length]
    in_file.close()
    return read_to_barcode_dict, total_reads

try:
    # New in Python 3.4
    from statistics import mean
except ImportError:
    def mean(list_of_values):
        """Calculate the mean average of a list of numbers."""
        # Quick and dirty, assumes already a list not an interator
        # so don't have to worry about getting the divisor.
        # Explicit float(...) to allow for Python 2 division.
        return sum(list_of_values) / float(len(list_of_values))

assert mean([1,2,3,4,5]) == 3

def get_fasta_stats(fasta):
    """function to get stats on a given fasta file"""
    with open(fasta, 'r') as seq:
        sizes = [len(record) for record in SeqIO.parse(seq, 'fasta')]
    min_contig = min(sizes)
    max_contig = max(sizes)
    avg_contig = mean(sizes)
    num_contig = len(sizes)
    return min_contig, max_contig, avg_contig, num_contig


def get_barcode_seq (read, in_fastq):
    """function to get the barcode seq from the fq file.
    fq file needs to be unzipped.
    currently assuming barcode if 8bp and at the start of each seq??
    """
    #THIS FUNCTION IS NOT USED.
    return True

def parse_line(line):
    """finction to parse line"""
    if not line.strip():
        return False #if the last line is blank
    if line.startswith("#"):
        return False
    #if not "Phytophthora" or "P." in line:
        #continue
    line = line.rstrip()
    out_put_str = ""
    cluster_line = line.rstrip("\n").split("\t")
    # are the database Phy singletons? - then we wont be interested
    number_of_reads_hitting_species = len(cluster_line)
    #set up some blank variables
    species = ""
    barcode_left = ""
    barcode_right = ""
    #database phy counter
    phy_count = 0
    return cluster_line, number_of_reads_hitting_species, species, barcode_left, \
           barcode_right, phy_count


def parse_tab_file_get_clusters(fasta, filename1, left, right, barcode_length, \
                                show_me_the_reads,\
                                read_prefix, out_file):
    """#script to open up a tab separeted clustering output and identify the
    species in the clustering"""
    #print coded_name_to_species_dict
    #call the functions
    left_read_to_barcode_dic, left_total_reads = turn_fq_to_dic(left, barcode_length)
    right_read_to_barcode_dic, right_total_reads = turn_fq_to_dic(right, barcode_length)
    
    #call the function to get fasta stats
    min_contig, max_contig, avg_contig, num_contig = get_fasta_stats(fasta)
    read_prefix = str(read_prefix)

    cluster_file = open (filename1, "r")
    summary_out_file = open(out_file, "w")
    #set so we dont replica the barcode output
    left_barcode_seen_set = set([])
    right_barcode_seen_set = set([])

    ITS_hitting_phy_count = 0
    #title for the results file
    if show_me_the_reads:
        title = "#species\tnumber_of_reads_hitting_species\tunique_barcodeR1\tunique_barcodeR2\treads_that_hit_species\n"
    else:
        title = "#species\tnumber_of_reads_hitting_species\tbarcodeR1\tbarcodeR2\n"
    summary_out_file.write(title)
    count = int(0)
    for line in cluster_file:
        count +=1
        #call function to parse line
        if not parse_line(line):
            continue
        reads_of_interest = ""
        cluster_line, number_of_reads_hitting_species, \
                      species, barcode_left, barcode_right, phy_count = parse_line(line)

        for member in cluster_line:
            member = member.rstrip()
            #is a memeber of the database in this cluster?
            if "Phytophthora" in member or "P." in member or "VHS" in member:
                # yes we are interested in this cluster
                phy_count = phy_count+1
                #print member
                species = species+member+" "
                #remove all the memebr which are database memebers
                number_of_reads_hitting_species = number_of_reads_hitting_species-1
        if number_of_reads_hitting_species >=1 and phy_count >0: #after the line, are there more memeber that are not database members?
            for member in cluster_line:
                if "Phytophthora" in member or "P." in member or "VHS" in member:
                    #no barcode for these - obviously!
                    continue
                if not member.startswith(read_prefix):
                    continue
                #print ("YYYEEEESSSSSS")
                #print cluster_line
                # we do not add the reads that hit the 
                reads_of_interest = reads_of_interest+member+" "
                #call func to get the bar codes
                #print "I am getting the left bar code for memebr ", member, " = ", left_read_to_barcode_dic[member]
                try:
                    #print "left_read_to_barcode_dic =",left_read_to_barcode_dic
                    left_barcode = left_read_to_barcode_dic[member]
                except:
                    left_barcode = "Not_found"
                # add to str
                if left_barcode not in left_barcode_seen_set:
                    left_barcode_seen_set.add(left_barcode)
                    barcode_left = barcode_left+left_barcode+" "
                #print "I am getting the right bar code for memebr", member, " = ", right_read_to_barcode_dic[member]
                try:
                    right_read_to_barcode_dic[member]
                    right_barcode = right_read_to_barcode_dic[member]
                except:
                    right_barcode = "Not_found"
                #add to str if we have not seen it yet. 
                if right_barcode not in right_barcode_seen_set:
                    right_barcode_seen_set.add(right_barcode)
                    barcode_right = barcode_right+right_barcode+" "

            ITS_hitting_phy_count = ITS_hitting_phy_count+number_of_reads_hitting_species
            #format the data
            if show_me_the_reads:
                data_output = "%s\t%d\t%s\t%s\t%s\n" %(species.rstrip(), number_of_reads_hitting_species,\
                                                       barcode_left,\
                                                       barcode_right, reads_of_interest)
            else:
                data_output = "%s\t%d\t%s\t%s\t\n" %(species.rstrip(), number_of_reads_hitting_species,\
                                                       barcode_left,\
                                                       barcode_right)
            #write out the data
            summary_out_file.write(data_output)
            phyto_read = (float(ITS_hitting_phy_count)/ num_contig)*100

    fasta_file_summary = """#Fasta file assembly summary: #min_contig = %d max_contig = %d avg_contig = %d
#Total number of assemblerd sequences = %d
#number of reads clustering with Phy = %d
#number of starting reads = %d
#percent of reads clustering with Phyto = %.2f""" %(min_contig, \
                                            max_contig, avg_contig, num_contig, ITS_hitting_phy_count,\
                                            right_total_reads, phyto_read)
    summary_out_file.write(fasta_file_summary)
    cluster_file.close()
    summary_out_file.close()
    return True

#############################################################################

#to run the script       

usage = """usage :

this scripts to summarise what clusters with what Phy species.


python parse_clusters.py -i clustering_outfile_already_decoded_from_temp_names -o summarise_clusters.out

command line option


"""

parser = OptionParser(usage=usage)

parser.add_option("-f","--fasta", dest="fasta", default=None,
                  help="assembled fasta file to get stats on")

parser.add_option("-i","--in", dest="in_file", default=None,
                  help="clustering out file")

parser.add_option("-l", "--left", dest="left", default="temp_not_trimmedr1.fq",
                  help="left reads unzipped fq file. Need this to get the barcode",
                  metavar="FILE")
parser.add_option("-r", "--right", dest="right", default="temp_not_trimmedr2.fq",
                  help="right reads unzipped fq file. Need this to get the barcode",
                  metavar="FILE")

parser.add_option("-b", "--barcode_length", dest="barcode_length", default=6,
                  help="length of the barcode used. Default 6 ",
                  metavar="FILE")

parser.add_option("--read_prefix", dest="read_prefix", default=None,
                  help="read_prefix from the fq file. Only needs "
                  " the first few letter e.g. read_prefix  M01157 ",
                  metavar="FILE")

parser.add_option("-s", "--show_me_the_reads", dest="show_me_the_reads", default=False,
                  help="show_me_the_reads in the output file for those that hit the Phy species."
                  " by default this is off, as the file could get very large. -s True if you want ...",
                  metavar="FILE")

parser.add_option("-o", "--out_prefix", dest="out_file", default="summarise_clusters.out",
                  help="prefix to the output filenames")


(options, args) = parser.parse_args()
#-f
fasta = options.fasta
#-i
in_file = options.in_file
# -l
left = options.left
# -r
right = options.right
# -b
barcode_length = options.barcode_length
# -o
out_file = options.out_file
# -s
show_me_the_reads = options.show_me_the_reads
#read_prefix
read_prefix = options.read_prefix

#run the program
barcode_length = int(barcode_length)

parse_tab_file_get_clusters(fasta, in_file, left, right, barcode_length, show_me_the_reads, \
                            read_prefix, out_file)

print "done"
