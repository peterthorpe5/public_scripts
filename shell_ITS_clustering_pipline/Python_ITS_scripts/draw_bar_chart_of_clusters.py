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
from math import log
#imports for graphs
#import seaborn as sns

import matplotlib
# this code added to prevent this error:
# self.tk = _tkinter.create(screenName, baseName, className, interactive, wantobjects, useTk, sync, use)
#_tkinter.TclError: no display name and no $DISPLAY environment variable
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import pylab
import numpy as np

# Turn off warning messages
import warnings
warnings.filterwarnings('ignore')


##############################################################################################
# drawing a histogram
#n, bins, patches = hist(data)

def parse_line(line):
    """function to parse a given line and return
    tab separated elements"""
    if not line.strip():
        return False #if the last line is blank
    if line.startswith("#"):
        return False
    cluster_line_split = line.rstrip("\n").split("\t")
    #print ("I assume the element are tab separated")
        #cluster_line_split = line.rstrip("\n").split()
    return cluster_line_split
    

def count_element_in_cluster(cluster_line_split):
    """func to count to total memebers in each cluster.
    Returns a dic with """
    species_set = set([])
    #counters to keep track
    members_count = 0
    species_count = 0
    for member in cluster_line_split:
        members_count = members_count+1
        try:
            if member.split("_")[1] =="sp":
                #split the name up to get the species name
                #this will be messed up by database enteries
               species_member = member.split("_")[2]
            else:    
                species_member = member.split("_")[1]
        except:
            ValueError
            species_member = member
        if member not in species_set:
            species_count = species_count+1
            species_set.add(species_member)
    return members_count, species_count

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

def covert_dict_to_list_of_value(in_dict):
    """function to convert a given dict, to convert it to a
    list of value. In dict should be something like this:
    count_dict= {1: 4, 2: 1}  """
    #need this value to plot the number of histogram bars
    number_of_keys = len(in_dict.keys())
    output_list = []
    key_list = []
    max_val = 0
    vals_for_bar_chart = []
    for key, val in in_dict.items():
        key_list.append(key)
	for i in range(0,val):
            output_list.append(key)
            #get the maximum val for graph
            if val > max_val:
                max_val = val
    for i in range(1, max(key_list)+1):
        try:
            val = in_dict[i]
            vals_for_bar_chart.append(val)
        except:
            KeyError
            vals_for_bar_chart.append(0)

    return sorted(output_list), number_of_keys, max_val, vals_for_bar_chart

def plot_hitstogram_graph(data_values, title, number_of_keys, max_val, file_in):
    """function to draw a histogram of a given list of values.
    http://matplotlib.org/1.3.0/examples/pylab_examples/histogram_demo_extended.html
    https://github.com/widdowquinn/Teaching-Data-Visualisation/blob/master/exercises/one_variable_continuous/one_variable_continuous.ipynb
    """
    
    #bins = max(data_values)
    #pylab.hist(data_values, facecolor='blue')
    pylab.hist(data_values, facecolor='green', alpha=0.6)
    pylab.grid(True)
    pylab.title(title+"_histogram")
    pylab.xlabel('number in cluster')
    pylab.ylabel('Count')
    pylab.savefig(file_in+"_"+title+'_histogram.png')
    plt.close()
    pylab.close()

    os.chdir('..')

def plot_individual_bar_chart_graph(data_values, title, number_of_keys, max_val, vals_for_bar_chart, file_in):
    """function to draw a bar of a given list of values.
    FOR these data this IS the correct type of graph.
    http://matplotlib.org/examples/api/barchart_demo.html
    https://github.com/widdowquinn/Teaching-Data-Visualisation/blob/master/exercises/one_variable_continuous/one_variable_continuous.ipynb
    bar(left, height, width=0.8, bottom=None, hold=None, **kwargs)
    """

    n_groups = len(vals_for_bar_chart)
    fig, ax = plt.subplots()
    index = np.arange(n_groups)
    bar_width = 0.9
    opacity = 0.4
    #print vals_for_bar_chart
    rects1 = plt.bar(index, vals_for_bar_chart, bar_width,
                     alpha=opacity, color='b') # label='whatever'
    plt.xlabel('number in cluster')
    plt.ylabel('Count')
    plt.title(title+"_barchart")
    plt.legend()
    pylab.grid(True)
    ax.set_yscale('symlog')
    ax.set_xscale('symlog')
    plt.tight_layout()
    plt.show()
    pylab.savefig(file_in+"_"+title+'_barchart.png')
    plt.close()
    pylab.close()


def plot_multi_bar_chart_graph(title1, vals_for_bar_chart1,\
                               title2, vals_for_bar_chart2,\
                               title3, vals_for_bar_chart3,\
                                    file_in):
    """function to draw a bar of a given list of values.
    FOR these data this IS the correct type of graph.
    http://matplotlib.org/examples/api/barchart_demo.html
    https://github.com/widdowquinn/Teaching-Data-Visualisation/blob/master/exercises/one_variable_continuous/one_variable_continuous.ipynb
    bar(left, height, width=0.8, bottom=None, hold=None, **kwargs)
    """
    n1_groups = len(vals_for_bar_chart1)
    n2_groups = len(vals_for_bar_chart2)
    n3_groups = len(vals_for_bar_chart3)
    fig = plt.figure(figsize=(10, 8), dpi=1200)
    ##    # Create subplot axes
    ax1 = fig.add_subplot(1, 3, 1)  # 1x3 grid, position 1
    ax2 = fig.add_subplot(1, 3, 2)  # 1x3 grid, position 2
    ax3 = fig.add_subplot(1, 3, 3)  # 1x3 grid, position 3
    
    index1 = np.arange(n1_groups)
    index2 = np.arange(n2_groups)
    index3 = np.arange(n3_groups)
    #print index
    bar_width = 0.9
    opacity = 0.6

    #graph1  
    rects1 = ax1.bar(index1, vals_for_bar_chart1, bar_width,
                     alpha=opacity,color='green') # label='whatever'
    ax1.set_xlabel('log number in cluster')
    ax1.set_ylabel('log Count')
    ax1.set_yscale('symlog')
    ax1.set_xscale('symlog')
    ax1.grid(True)
    ax1.set_title(title1+"_barchart")
    
    # graph 2
    rects2 = ax2.bar(index2, vals_for_bar_chart2, bar_width,
                     alpha=opacity,color='blue') # label='whatever'
    ax2.set_xlabel('log number in cluster')
    ax2.set_ylabel('log Count')
    ax2.set_yscale('symlog')
    ax2.set_xscale('symlog')
    ax2.grid(True)
    ax2.set_title(title2+"_barchart")
   
    # graph 3
    rects3 = ax3.bar(index3, vals_for_bar_chart3, bar_width,
                     alpha=opacity,color='pink') # label='whatever'
    ax3.set_xlabel('test cluster')
    ax3.set_ylabel('test_count_Count')
    pylab.grid(True)
    ax3.set_title(title3+"_barchart")
    fig.tight_layout()
    fig
    pylab.savefig(file_in+'_GROUP_barchart.png')
    pylab.close()


    
def parse_tab_file_get_clusters(filename1, out_file):
    """#script to open up a tab separeted clustering output and identify the
    species in the clustering"""
    cluster_file = open (filename1, "r")
    summary_out_file = open(out_file, "w")
    #dictionaries for keeping the counts
    member_in_cluster_to_count_dict = dict()
    species_in_cluster_count_dict = dict()

    # a way of keeping track of the iteration
    interation_count = int(0)
    #iterate through the file
    for line in cluster_file:
        interation_count +=1
        #call the func to split up the line
        cluster_line_split = parse_line(line.rstrip())
        if not cluster_line_split:
            #this could be a blank line or a line that starts with #
            continue
        
        # call the function to get the number of elements and spceis. 
        members_count, species_count = count_element_in_cluster(cluster_line_split)
        try:
            #if we have seen this count before, then just add one to it.
            member_in_cluster_to_count_dict[members_count] +=1
        except:
            KeyError
            #print ("not seen this before")
            #not seen this before, set up a new dic element and make the equal 1
            member_in_cluster_to_count_dict[members_count] = 1
        try:
            #if we have seen this count of species before,
            #then just add one to it.
            species_in_cluster_count_dict[species_count] +=1
        except:
            KeyError
            species_in_cluster_count_dict[species_count] = 1


    #call the function to convert dic to list
    #print "species_in_cluster_count_dict: ", species_in_cluster_count_dict
    species_in_cluster_list, species_number_of_keys, species_max_val, \
                             species_vals_for_bar_chart = covert_dict_to_list_of_value(species_in_cluster_count_dict)

    #print "member_in_cluster_to_count_dict: ", member_in_cluster_to_count_dict
    member_in_cluster_list, member_number_of_keys, member_max_val, \
                            member_vals_for_bar_chart = covert_dict_to_list_of_value(member_in_cluster_to_count_dict)

    #plot one graph per figure
    #plot_individual_bar_chart_graph(species_in_cluster_list, "species_in_cluster",\
                         #species_number_of_keys, species_max_val, \
                         #species_vals_for_bar_chart, filename1)
    
    #plot_individual_bar_chart_graph(member_in_cluster_list, "total_member_in_cluster", \
                         #member_number_of_keys, member_max_val, \
                         #member_vals_for_bar_chart, filename1)



    #plot histogram - inappropriate for the data
    #plot_hitstogram_graph(species_in_cluster_list, "species_in_cluster",\
                         #species_number_of_keys, species_max_val, filename1)
    
    #plot_hitstogram_graph(member_in_cluster_list, "total_members_in_cluster", \
                         #member_number_of_keys, member_max_val, filename1)
    #plot_multi_bar_chart_graph
    plot_multi_bar_chart_graph("species in cluster",species_vals_for_bar_chart,\
                               "total member in cluster", member_vals_for_bar_chart,\
                               "test", member_vals_for_bar_chart,\
                               filename1)


    return True

#############################################################################

#to run the script       

usage = """usage :

script to graphically represent clustering output.


python draws_bar...py -i clustering_file -o summarise_clusters.out

requires:
use pip install ...
 seaborn
 matplotlib
 numpy

"""

parser = OptionParser(usage=usage)

parser.add_option("-i","--in", dest="in_file", default=None,
                  help="clustering out file")
parser.add_option("--heatmap", dest="heatmap", default=False,
                  help="draw a heat map of the species clustering")

parser.add_option("-o", "--out_prefix", dest="out_file", default="summarise_clusters.out",
                  help="prefix to the output filenames")


(options, args) = parser.parse_args()

in_file = options.in_file
heatmap = options.heatmap
out_file = options.out_file

#run the program

parse_tab_file_get_clusters(in_file, out_file)

if heatmap:
    from draw_heat_map import *

print "done"
