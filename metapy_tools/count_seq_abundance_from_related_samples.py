#!/usr/bin/env python
# title: duplication_number_to_filename
# (c) The James Hutton Institute 2018
# Author: Peter Thorpe, Leighton Pritchard
import sys
import os
import argparse
import subprocess
import operator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
# This line will not work for me on my set up.
#import matplotlib.pyplot as plt
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
from collections import OrderedDict
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

python count_seq_abundance_from_related_samples.py

this will open in two assembled fasta files, from two comparable
experiments and count the amount of times THAT seq is seen.

require: Biopython, matplotlib, pandas

"""
if "--help" or "-h" in sys.argv:
    print(usage)


def get_args():
    parser = argparse.ArgumentParser(description=" count_seq_abundance_from_related_samples.py",
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("count_seq_abundance_from_related_samples")[0]
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument("-s", "--suffix", dest='suffix',
                          action="store",
                          default=".assembled.fastq.bio.chopped.fasta",
                          type=str,
                          help="suffix: .assembled.fastq.bio.chopped.fasta")
    optional.add_argument("-1", "--in1", dest='in1',
                          action="store",
                          default=None,
                          type=str,
                          help="infile1")

    optional.add_argument("-2", "--in2", dest='in2',
                          action="store",
                          default=None,
                          type=str,
                          help="infile2")

    optional.add_argument("-o", "--out", dest='out',
                          action="store",
                          default="results.txt",
                          type=str,
                          help="outfile name")

    optional.add_argument('--version',
                          action='version',
                          version="%s: abundance.py " + VERSION)
    args = parser.parse_args()
    return args, file_directory


def run_parse_cmd(cmd):
    """func to run the parser command as this may be called a few times"""
    # logger.info("%s make cmd_parse command", cmd)
    #  removed check-True
    #print(cmd)
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)


def math_log(invalue):
    """try to log the value, if 0 or 1 return the value"""
    if invalue == 0:
        return invalue
    return log(invalue, 10)

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
                if x.startswith(args.in1):
                    wanted_file1 = (os.path.join(dirpath, x))
                    print(wanted_file1)
                    continue
    for dirpath, subdirs, files in os.walk(directory):
        for j in files:
            if j.endswith(SUFFIX):
                if j.startswith(args.in2):
                    wanted_file2 = (os.path.join(dirpath, j))
                    print(wanted_file2)
                    both_files = "%s_%s.temp.fasta" % (os.path.split(wanted_file1)[-1].split("_L")[0],
                                                       j.split("_L")[0])
                    cmd = " ".join(["cat",
                                    wanted_file1,
                                    wanted_file2,
                                    ">",
                                    both_files])
                    run_parse_cmd(cmd)
                    continue
    infile_1_count = defaultdict(int)
    infile_2_count = defaultdict(int)
    for seq_record in SeqIO.parse(both_files, "fasta"):
        infile_1_count[str(seq_record.seq)] = 0
        infile_2_count[str(seq_record.seq)] = 0

    for seq_record in SeqIO.parse(wanted_file1, "fasta"):
        infile_1_count[str(seq_record.seq)] += 1

    for seq_record in SeqIO.parse(wanted_file2, "fasta"):
        infile_2_count[str(seq_record.seq)] += 1

    infile1_seq_count = []
    infile2_seq_count = []
    sequneces = []
    interative_count = 0
    #print("infile_1_count = ", infile_1_count)
    # sort the dictionary by its values (i think?)
    #sorted_infile_1_count = OrderedDict(sorted(d.items(), key=lambda x: x[1]))
    for sequence, count in sorted(infile_1_count.items(), key=lambda x: x[1], reverse=True):
        interative_count = interative_count + 1
        sequneces.append(sequence)
        #infile1_seq_count.append(math_log(count))
        infile1_seq_count.append(count)
        count2 = infile_2_count[sequence]
        #infile2_seq_count.append(math_log(count2))
        infile2_seq_count.append(count2)

    # print("length = infile1_seq_count", len(infile1_seq_count))
    # print("length = infile2_seq_count", len(infile1_seq_count))
    # print("count = ", interative_count)

    # example used: https://python-graph-gallery.com/124-spaghetti-plot/
    #df=pd.DataFrame({'y': range(0,interative_count),'x1': infile1_seq_count, 'x2':infile2_seq_count})
    df=pd.DataFrame({'x': range(0,interative_count), 'y1': infile1_seq_count, 'y2': infile2_seq_count})
    print(df)
    plt.style.use('seaborn-darkgrid')
    # multiple line plot
    num=0
    for column in df.drop('x', axis=1):
        num+=1
        plt.plot(df['x'], df[column], marker='',
                 linewidth=0.3, alpha=0.9, label=column)

        # Add legend
    plt.legend(loc=2, ncol=2)

    # Add titles
    name = args.in1.split("_L")[0] + "_and_" + args.in2.split("_L")[0]
    plt.title("Graph representing the number of times an assembled sequence is seen for %s" % name,
              fontsize=10)
    plt.xlabel("LOG Number of of times a sequences occurs")
    plt.ylabel("Number of sequences")
    plt.tick_params(axis='both', which='major', labelsize=8)
    plt.tick_params(axis='both', which='minor', labelsize=4)
    plt.xticks(fontsize=4, rotation=90)
    pylab.grid(True)
    plt.show()
    plt.savefig(name + "_sequences_seen_plot" + ".png")




