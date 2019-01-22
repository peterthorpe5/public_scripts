#!/usr/bin/env python2.7
# title: get blast stats and draw a graph of the results

# author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

# imports

import sys
import os
from optparse import OptionParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
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


def read_tab_file(tab_output):
    """read in the tab file. Reads whole file into memory.
    Could be altered for more efficiency
    """
    with open(tab_output) as file:
        return file.read().split("\n")


def parse_line(line):
    """function to parse a given line and return
    tab separated elements"""
    if not line.strip():
        return False #if the last line is blank
    if line.startswith("#"):
        return False
    line_split = line.rstrip("\n").split("\t")
    #print ("I assume the element are tab separated")
        #cluster_line_split = line.rstrip("\n").split()
    return line_split


def get_perc(data, perc):
    """func to return a list of perct identity from
    coloumn 3 in the tab blast file.
    Take in a tab sep blast line and the list which is
    appends to"""
    perc.append(int(round(float(data[2]))))
    return  perc


def get_bit_list(data, bit):
    """func to return a list of bit scores from
    coloumn 3 in the tab blast file.
    Take in a tab sep blast line and the list which is
    appends to"""
    bit.append(int(round(float(data[11]))))
    return bit


def get_alignmemt_list(data, align):
    """func to return a list of align length from
    coloumn 3 in the tab blast file.
    Take in a tab sep blast line and the list which is
    appends to"""
    align.append(int(round(float(data[3]))))
    return align


def plot_hitstogram_graph(data_values, title, file_in):
    """function to draw a histogram of a given list of values.
    http://matplotlib.org/1.3.0/examples/pylab_examples/
    histogram_demo_extended.html
    https://github.com/widdowquinn/Teaching-Data-Visualisation/
    blob/master/exercises/one_variable_continuous/
    one_variable_continuous.ipynb
    """

    #bins = max(data_values)
    #pylab.hist(data_values, facecolor='blue')
    pylab.hist(data_values, facecolor='green', alpha=0.6)
    pylab.grid(True)
    pylab.title(title + "_histogram")
    pylab.xlabel('Percentage Identity')
    pylab.ylabel('Number in Bin')
    #pylab.savefig(file_in + "_" + title + '_histogram.png')
    pylab.savefig(file_in + "_" + title + '_histogram.pdf', format=pdf)
    plt.close()
    pylab.close()
    os.chdir('.')


def plot_multi_histogram_graph(title1, vals_for_hist1,
                               title2, vals_for_hist2,
                               title3, vals_for_hist3,
                               file_in):
    """function to draw a histo of a given list of values.
    FOR these data this IS the correct type of graph.
    http://matplotlib.org/examples/api/barchart_demo.html
    https://github.com/widdowquinn/
    Teaching-Data-Visualisation/blob/master/
    exercises/one_variable_continuous/
    one_variable_continuous.ipynb
    bar(left, height, width=0.8, bottom=None,
    hold=None, **kwargs)
    """

    fig = plt.figure(figsize=(10, 8), dpi=1200)
    #    # Create subplot axes
    ax1 = fig.add_subplot(1, 3, 1)  # 1x3 grid, position 1
    ax2 = fig.add_subplot(1, 3, 2)  # 1x3 grid, position 2
    ax3 = fig.add_subplot(1, 3, 3)  # 1x3 grid, position 3
    matplotlib.rcParams.update({'font.size': 8})

    # print (index)
    bar_width = 0.9
    opacity = 0.6

    # graph1 pylab.hist
    rects1 = ax1.hist(vals_for_hist1,
                      facecolor='green',
                      alpha=0.6) # label='whatever'
    ax1.set_xlabel('Percentage Identity')
    ax1.set_ylabel('Number in Bin')
    #ax1.set_yscale()
    #ax1.set_xscale()
    ax1.grid(True)
    ax1.set_title(title1)

    # graph 2
    rects2 = ax2.hist(vals_for_hist2,
                      facecolor='blue',
                      alpha=0.6)  # label='whatever'
    ax2.set_xlabel('Bit Score')
    ax2.set_ylabel('Number in Bin')
    #ax2.set_yscale()
    #ax2.set_xscale()
    ax2.grid(True)
    ax2.set_title(title2)

    # graph 3
    rects3 = ax3.hist(vals_for_hist3,
                      facecolor='red',
                      alpha=0.6)  # label='whatever'
    ax3.set_xlabel('Alignmnet Length')
    ax3.set_ylabel('Number in Bin')
    pylab.grid(True)
    ax3.set_title(title3 + "_histogram")
    fig.tight_layout()
    fig
    pylab.savefig(file_in + '_histogram.pdf')
    pylab.close()


def parse_tab_file(tab_output, out_file):
    """#script to open up a tab blast output and plot the
    percentage identity."""
    blast = read_tab_file(tab_output)
    perc = []
    bit = []
    align = []
    for line in blast:
        if parse_line(line) == False:
            continue
        else:
            data = parse_line(line)
            perc = get_perc(data, perc)
            bit = get_bit_list(data, bit)
            align = get_alignmemt_list(data, align)
    # plot_hitstogram_graph(perc, "percentage_identity",
                          # "gene_models_test")

    plot_multi_histogram_graph("percentage_Identity",
                               perc,
                               "Bits_Scores",
                               bit,
                               "Alignment_lengths",
                               align,
                               out_file)

    summary_out_file = open(out_file, "w")



if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:
$ python blast_stats.py -i blast.tab -o outfile
script to plot blast tab results.
"""

parser = OptionParser(usage=usage)
parser.add_option("-i", dest="in_file",
                  default=None,
                  help="blast tabular output")

parser.add_option("-o", "--out", dest="out",
                  default=None,
                  help="Output filename",
                  metavar="FILE")

(options, args) = parser.parse_args()

in_file = options.in_file
out = options.out


# Run as script
if __name__ == '__main__':
    if not os.path.isfile(in_file):
        print("sorry cannot find you %s file" % in_file)
        os._exit(0)  
    parse_tab_file(in_file, out)

