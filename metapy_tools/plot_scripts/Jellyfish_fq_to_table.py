#!/usr/bin/env python
# title: duplication_number_to_filename
# (c) The James Hutton Institute 2018
# Author: Peter Thorpe, Leighton Pritchard
import sys
import os
import argparse
import subprocess
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

python Jellyfish_fq_to_table.py

this will output a file for each file in the sub directories that match
*.assembled.fastq , as ouputted by metapy. By default these are
deleted by metapy so you have to run with --cleanup no


require: matplotlib, and jellyfish to be in your path

"""
if "--help" or "-v" in sys.argv:
    print(usage)


def get_args():
    parser = argparse.ArgumentParser(description="jellyfish kmer count ",
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("Jelly")[0]
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument("-s", "--suffix", dest='suffix',
                          action="store",
                          default=".assembled.fastq.bio.chopped",
                          type=str,
                          help="suffix: .assembled.fastq.bio.chopped")
    
    optional.add_argument("-k", "--kmer", dest='kmer',
                          action="store",
                          default="11",
                          type=str,
                          help="kmer size for jellyfish to count")

    optional.add_argument("-t", "--threads", dest='threads',
                          action="store",
                          default="4",
                          type=str,
                          help="kmer size for jellyfish to count")

    optional.add_argument("-o", "--out", dest='out',
                          action="store",
                          default="results.txt",
                          type=str,
                          help="outfile name")

    optional.add_argument("-v", '--version',
                          action='version',
                          version="%s: Jellyfish_fq.py " + VERSION)
    args = parser.parse_args()
    return args, file_directory


def run_parse_cmd(cmd):
    """func to run the parser command as this may be called a few times"""
    # logger.info("%s make cmd_parse command", cmd)
    #  removed check-True
    print(cmd)
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)

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
                #f_out = open(PREFIX + "_DEREP_NUMBERS.txt", "w")
                j_cmd = " ".join(["jellyfish",
                                  "count",
                                  "-C",
                                  "-m",
                                  "%s" % args.kmer,
                                  "-s",
                                  "10000000",
                                  "-t",
                                  "%s" % args.threads,
                                  "%s" % (os.path.join(dirpath, x)),
                                  "-o",
                                  "%s_reads.jf" % PREFIX])
                run_parse_cmd(j_cmd)
                h_cmd = " ".join(["jellyfish",
                                  "histo",
                                  "-t",
                                  "%s" % args.threads,
                                  "%s_reads.jf" % PREFIX,
                                  ">",
                                  "%s_reads_%s.histo" % (PREFIX, args.kmer)])
                run_parse_cmd(h_cmd)
                rm_cmd = " ".join(["rm",
                                   "%s_reads.jf" % PREFIX])
                run_parse_cmd(rm_cmd)
                
                                
                                   
##                derep_list = []
##                derep_position = []
##                counts = []
##                duplication_count_dict = defaultdict(int)
##                for seq_record in SeqIO.parse(os.path.join(dirpath, x), "fasta"):
##                    derep_nums = seq_record.id.split("_")[1] + "\n"
##                    duplication_count_dict[int(derep_nums.rstrip())] += 1
##                    f_out.write(derep_nums)
##                    #derep_list.append(int(derep_nums.rstrip()))
##                for derep_nums, count in sorted(duplication_count_dict.items()):
##                    derep_list.append(log(derep_nums, 10))
##                    counts.append(log(count, 10))                   
##                    
##                df=pd.DataFrame({'xvalues': counts, 'yvalues': derep_list})
##                plt.plot( 'xvalues', 'yvalues', data=df, linestyle='', marker='o',
##                          markersize=2)
##                plt.xlabel('LOG Number of unique sets of this size found')
##                plt.ylabel('LOG Number of reads in rereplicated read set')
##                plt.tick_params(axis='both', which='major', labelsize=8)
##                plt.tick_params(axis='both', which='minor', labelsize=4)
##                plt.xticks(fontsize=4, rotation=90)
####                #sns.kdeplot(df.x, df.y, cmap="Reds", shade=True)
####                pylab.hist(derep_list, bins=20, histtype="bar", facecolor='blue', alpha=0.6)
####                pylab.xlabel('NUmber of unique reads')
####                pylab.ylabel('Number in Bin')
##                pylab.grid(True)
##                plt.title("Unique reads for: %s" % PREFIX )
##                plt.show()
##                plt.savefig(PREFIX + ".png")
##                    # print derep_nums
##                f_out.close()
                    



