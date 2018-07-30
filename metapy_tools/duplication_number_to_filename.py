#!/usr/bin/env python
# title: duplication_number_to_filename
# (c) The James Hutton Institute 2018
# Author: Peter Thorpe, Leighton Pritchard
import sys
import os
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import matplotlib
# this code added to prevent this error:
# self.tk = _tkinter.create(screenName, baseName,
# className, interactive, wantobjects, useTk, sync, use)
# _tkinter.TclError: no display name and
# no $DISPLAY environment variable
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
#import seaborn as sns
plt.style.use('seaborn')
import pylab
from math import log



VERSION = "summerise results: v0.01"
if "--version" in sys.argv:
    print(VERSION)
    sys.exit(1)

usage = """

python -h

python duplication_number_to_filename.py

this will output a file for each file in the sub directories that match
*astqfor_swarm.fasta , as ouputted by metapy. By default these are
deleted by metapy so you have to run with --cleanup no

in the output file the number of dureplicated reads will be output.

e.g.

>c6533710b662119106c74f34c3374753_3

this was forms by the deduplication of a read set of which 3 were identical. 

require: Biopython, matplotlib, pandas

"""
if "--help" or "-h" in sys.argv:
    print(usage)


def get_args():
    parser = argparse.ArgumentParser(description="dereplication numbers ",
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("duplication_number_to_filename")[0]
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument("-s", "--suffix", dest='suffix',
                          action="store",
                          default="astqfor_swarm.fasta",
                          type=str,
                          help="suffix astqfor_swarm.fasta")

    optional.add_argument("-o", "--out", dest='out',
                          action="store",
                          default="results.txt",
                          type=str,
                          help="outfile name")

    optional.add_argument('--version',
                          action='version',
                          version="%s: metapy.py " + VERSION)
    args = parser.parse_args()
    return args, file_directory


args, FILE_DIRECTORY = get_args()
# Run as script
if __name__ == '__main__':
    SUFFIX = str(args.suffix)
    # call the function to get a list of results wanted
    directory = "."
    # get all folder names os.walk(directory)
    for dirpath, subdirs, files in os.walk(directory):
        for x in files:
            #print result
            if x.endswith(SUFFIX):
            #if os.path.split(result)[-1].endswith(SUFFIX):
                PREFIX = x.split(".assembled")[0]
                # print os.path.join(dirpath, x), SUFFIX
                f_out = open(PREFIX + "_DEREP_NUMBERS.txt", "w")
                derep_list = []
                derep_position = []
                counts = []
                duplication_count_dict = defaultdict(int)
                for seq_record in SeqIO.parse(os.path.join(dirpath, x), "fasta"):
                    derep_nums = seq_record.id.split("_")[1] + "\n"
                    duplication_count_dict[int(derep_nums.rstrip())] += 1
                    f_out.write(derep_nums)
                    #derep_list.append(int(derep_nums.rstrip()))
                for derep_nums, count in sorted(duplication_count_dict.items()):
                    derep_list.append(log(derep_nums, 10))
                    counts.append(log(count, 10))                   
                    
                df=pd.DataFrame({'xvalues': counts, 'yvalues': derep_list})
                plt.plot( 'xvalues', 'yvalues', data=df, linestyle='', marker='o',
                          markersize=2)
                plt.xlabel('LOG Number of unique sets of this size found')
                plt.ylabel('LOG Number of reads in rereplicated read set')
                plt.tick_params(axis='both', which='major', labelsize=8)
                plt.tick_params(axis='both', which='minor', labelsize=4)
                plt.xticks(fontsize=4, rotation=90)
##                #sns.kdeplot(df.x, df.y, cmap="Reds", shade=True)
##                pylab.hist(derep_list, bins=20, histtype="bar", facecolor='blue', alpha=0.6)
##                pylab.xlabel('NUmber of unique reads')
##                pylab.ylabel('Number in Bin')
                pylab.grid(True)
                plt.title("Unique reads for: %s" % PREFIX )
                plt.show()
                plt.savefig(PREFIX + ".png")
                    # print derep_nums
                f_out.close()
                    



