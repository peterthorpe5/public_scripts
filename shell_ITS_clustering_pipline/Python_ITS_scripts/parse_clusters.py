#!/usr/bin/env python
#title: Parse clusters and find the corresponding spcies
#in the cluster
#author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

#imports
import os
import sys
from optparse import OptionParser
import datetime
import os
from sys import stdin,argv


##############################################################################################
def coded_name_to_species(database_file):
    """functiong takes the already generated tab separated
    database of coded name to species file. Returns a dic
    of coded_name to species"""
    with open(database_file) as file:
        data= file.read().split("\n")
    coded_name_to_species_dict=dict()
    for line in data:
        if not line.strip():
            continue #if the last line is blank
        if line.startswith("#"):
            continue 
        #print line
        coded_name, species = line.split("\t")
        #print coded_name, species
        coded_name_to_species_dict[coded_name.rstrip()] = species.rstrip()
    return coded_name_to_species_dict

def rev_coded_name_to_species(database_file):
    """functiong takes the already generated tab separated
    database of coded name to species file. Returns a dic
    of coded_name to species"""
    with open(database_file) as file:
        data= file.read().split("\n")
    rev_coded_name_to_species_dict=dict()
    for line in data:
        if not line.strip():
            continue #if the last line is blank
        if line.startswith("#"):
            continue 
        #print line
        coded_name, species = line.split("\t")
        #print coded_name, species
        rev_coded_name_to_species_dict[species.rstrip()] = coded_name.rstrip()
    return rev_coded_name_to_species_dict

def parse_tab_file_get_clusters(filename1, database, out_file):
    """#script to open up a tab separeted clustering output and identify the
    species in the clustering"""
    #call the function to get the dictionary populated with the database
    coded_name_to_species_dict = coded_name_to_species(database)
    rev_coded_name_to_species_dic = rev_coded_name_to_species(database)
    #print coded_name_to_species_dict
    cluster_file = open (filename1, "r")
    summary_out_file = open(out_file, "w")

    count = int(0)
    for line in cluster_file:
        if not line.strip():
            continue #if the last line is blank
        if line.startswith("#"):
            continue
        out_put_str = ""
        if "\t" in line:
            cluster_line = line.rstrip("\n").split("\t")
        else:
            cluster_line = line.rstrip("\n").split()
        count +=1
        for member in cluster_line:
            #print member
            try:
                species = coded_name_to_species_dict[member]
                #print species
            except:
                KeyError
                species = member
            try:
                species = rev_coded_name_to_species_dic[member]
                #print species
            except:
                KeyError
                species = member
            #add the info to a str which we will write at the end of the cluster line    
            cluster_summary = "%s\t" %(species)
            out_put_str = out_put_str+cluster_summary
        summary_out_file.write(out_put_str+"\n")

    cluster_file.close()
    summary_out_file.close()
    return True

#############################################################################

#to run the script       

usage = """usage :

this scripts return a 'master file' summarising what is in the clusters.

What is in the clusters is determined by a pre-made databse file of "coded" seq name
to species ...

python parse_clusters.py -i clustering_outfile -d name_to_species_database.txt -o summarise_clusters.out

command line option

"""

parser = OptionParser(usage=usage)

parser.add_option("-i","--in", dest="in_file", default=None,
                  help="clustering out file")

parser.add_option("-d", "--database", dest="database", default="name_to_species_database.txt",
                  help="prefix to alter the id names",
                  metavar="FILE")
parser.add_option("-o", "--out_prefix", dest="out_file", default="summarise_clusters.out",
                  help="prefix to the output filenames")


(options, args) = parser.parse_args()

in_file = options.in_file
database = options.database
out_file = options.out_file

#run the program

parse_tab_file_get_clusters(in_file, database, out_file)

print "done"
