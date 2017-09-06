#!/usr/bin/env python
# Code to iterogate clustering with a list
#
# (c) The James Hutton Institute 2016-2017
# Author: Peter Thorpe

import os
from sys import stdin,argv
import sys
from optparse import OptionParser
from collections import Counter
import collections

if "-v" in sys.argv or "--version" in sys.argv:
    print("0.01 - get the clusters from a list of seq of interest")
    sys.exit(os.system(cmd))


def parse_clusters(clusters):
    """funct to return list of cluserts"""
    with open(clusters) as handle:
        return handle.read().split("\n")


def get_set_of_interest(infile):
    """funtcion to load in a list of gene names"""
    with open(infile) as handle:
        data = handle.read().split()
        outset = set([])
        for entry in data:
            outset.add(entry.rstrip())
        return outset


##################################################################
usage = """Use as follows:

$ python unique_comprisons.py -i list_of_gene -c cluster_file -a allfile.wanted -o outfile

"""

parser = OptionParser(usage=usage)

parser.add_option("-i","--wanted", dest="infile",
                  default=None,
                  help="infile with the names of interest",
                  metavar="FILE")

parser.add_option("-c","--clusters", dest="clusters",
                  default="Orthofinder_OrthologousGroups_final.txt",
                  help="clusters file",
                  metavar="FILE")
parser.add_option("-a","--all", dest="all_file",
                  default="all_unique_v1.0.txt",
                  help="all unique gene names file",
                  metavar="FILE")
parser.add_option("-o", "--out", dest="out", default="result.out",
                  help="output filenames")


(options, args) = parser.parse_args()

infile = options.infile
clusters = options.clusters
all_file = options.all_file
out = options.out
################################################################

if __name__ == '__main__':
    if not os.path.isfile(clusters):
        print("sorry cannot find you %s file" % clusters)
        os._exit(0)
    if not os.path.isfile(clusters):
        print("sorry cannot find you %s infile" % infile)
        os._exit(0)
    working_dir = os.getcwd()
    dest_dir = os.path.join(working_dir, 'results')
    try:
        os.makedirs(dest_dir)
    except OSError:
        print("folder already exists, I will write over what is in there!!")
    cluster_data = parse_clusters(clusters)
    wanted = get_set_of_interest(infile)
    all_unique = get_set_of_interest(all_file)

    # track the interesting clusters so we dont get repeats
    clusters_of_interest = set([])
    outfile_path = os.path.join(working_dir, 'results', out)
    f_out = open(outfile_path, "w")
    allowed = ['Mpe', 'Mca', 'Api','Rpa', 'Dno']
    print "starting wanted list = %d " % len(wanted)

    # parser through the cluster file
    wanted_tracked_count = 0
    total_unique_counter = Counter({'Mpe':0, 'Mca':0, 'Api': 0,
                                   'Rpa':0, 'Dno':0})
    Total_elements_counter = Counter({'Mpe':0, 'Mca':0, 'Api': 0,
                                   'Rpa':0, 'Dno':0})
    Total_elements_matching_wanted_counter = Counter({'Mpe':0, 'Mca':0, 'Api': 0,
                                                      'Rpa':0, 'Dno':0})
    for line in cluster_data:
        line = line.rstrip()
        unique_counter_counter = Counter({'Mpe':0, 'Mca':0, 'Api': 0,
                                         'Rpa':0, 'Dno':0})
        species_counter = Counter({'Mpe':0, 'Mca':0, 'Api': 0,
                                   'Rpa':0, 'Dno':0})
        cluster_elements = line.split()
        # each entry separetly
        for gene in cluster_elements:
            gene = gene.rstrip()
            # only if the cluster contains a wanted gene
            prefix = gene[:3]
            Total_elements_counter[prefix] += 1
            if gene in wanted:
                # check to see if we have seen this line before
                prefix = gene[:3]
                Total_elements_matching_wanted_counter[prefix] += 1
                if not line in clusters_of_interest:
                    clusters_of_interest.add(line.rstrip())
                    # count through th cluster again to see what speices are there
                    for gene in cluster_elements:
                        gene = gene.rstrip()
                        prefix = gene[:3]
                        if prefix in allowed:
                            # double check only allowed species are counted
                            species_counter[prefix] += 1
                        if gene in all_unique:
                            unique_counter_counter[prefix] += 1
                            total_unique_counter[prefix] += 1
                            if gene in wanted:
                                wanted_tracked_count = wanted_tracked_count + 1
                    #print len(line.split())
                    #print species_counter
                    #print unique_counter_counter
                    extra = "Cluster size = \t"
                    species_counter_od = collections.OrderedDict(sorted(species_counter.items()))
                    species_counter_od = collections.OrderedDict(sorted(unique_counter_counter.items()))

                    out_formatted = "%s%d\t\tSPECIES: %s\t\tUNIQUE:\t%s\n" % (extra,
                                                                          len(line.split()),
                                                                          species_counter_od,
                                                                          species_counter_od)
                    f_out.write(out_formatted)
    print "total found = %d" % wanted_tracked_count
    print "total_unique_counter = ", collections.OrderedDict(sorted(total_unique_counter.items()))
    print "Total_elements_counter = ", collections.OrderedDict(sorted(Total_elements_counter.items()))
    print "Total_elements_matching_wanted_counter = ", collections.OrderedDict(sorted(Total_elements_matching_wanted_counter.items()))
    f_out.close()
    f_out.close()

