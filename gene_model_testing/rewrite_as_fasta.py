#script to re-write badly formatted fasta file. Remove duplicates,
#or get seq or interest.

import sys
import os
import re
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
    return sizes, min_contig, max_contig, avg_contig, num_contig


def reformat_as_fasta(outfile):
    "this function re-write a file as a fasta file"
    f= open(outfile, 'w')
    for seq_record in SeqIO.parse("temp_fa.fa", "fasta"):
        # replace stop codons as diamond does not like these
        seq = str(seq_record.seq)
        # seq = seq.replace("*", "")
        seq_record.seq = Seq(re.sub('[^A-Za-z]',"",str(seq).upper()))
        seq_record.description=""
        SeqIO.write(seq_record, f, "fasta")
    f.close()
    return True


def plot_multi_histogram_graph(vals, file_in):
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

    matplotlib.rcParams.update({'font.size': 6})

    # Create subplot axes
    ax1 = fig.add_subplot(1, 4, 1)  # 1x4 grid, position 1
    ax2 = fig.add_subplot(1, 4, 2)  # 1x4 grid, position 2
    ax3 = fig.add_subplot(1, 4, 3)
    ax4 = fig.add_subplot(1, 4, 4)

    bar_width = 0.9
    opacity = 0.6

    # graph1 pylab.hist
    rects1 = ax1.hist(vals,
                      facecolor='green',
                      alpha=0.6) # label='whatever'
    ax1.set_ylabel('Gene Sizes')
    # ax1.set_yscale()
    # ax1.set_xscale()
    ax1.grid(True)
    ax1.set_title("Histogram of Gene Sizes")

    # graph2  pylab.hist
    rects1 = ax2.hist(vals,
                      facecolor='green',
                      alpha=0.6) # label='whatever'
    ax2.set_ylabel('Gene Sizes')
    ax2.set_yscale("log")
    # ax1.set_xscale()
    ax2.grid(True)
    ax2.set_title("Histogram of Log Gene Sizes")


    # graph 3
    rects2 = ax3.boxplot(vals, 0, 'gD')
    # ax3.set_xlabel()
    ax3.set_ylabel('Gene Sizes')
    # ax2.set_yscale()
    # ax2.set_xscale()
    ax3.grid(True)
    ax3.set_title("Boxplot of Gene Sizes")

    # graph 4
    rects2 = ax4.boxplot(vals, 0, 'gD')
    # ax4.set_xlabel()
    ax4.set_ylabel('Log Gene Sizes')
    ax4.set_yscale("log")
    # ax2.set_xscale()
    ax4.grid(True)
    ax4.set_title("Boxplot of Log Gene Sizes")
    fig.tight_layout()
    fig
    pylab.savefig(file_in.split(".fa")[0] + '.pdf')
    pylab.close()


def reformat_as_force(filename, length):
    """this function re-write a file as a fasta file.
    This used to use Biopython, but something kept breaking
    it. So, using a nasty parser instead."""
    f = open("temp_fa.fa", 'w')
    with open(filename) as handle:
        seq = ""
        count = 0
        name_set = set([])
        for line in handle:
            if line.startswith(">"):
                if count > 0:
                    if len(seq.replace(".", "")) > int(length):
                        data = "%s%s" % (name, seq.replace(".", ""))
                        f.write(data)
                    seq = ""
                name = line
                count = count + 1
            else:
                seq = seq + line
    f.close()


if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

converts

$ python rewrite_as_fasta.py -i in.fasta -l min_length_of_seq (default(3)) --not_wanted --wanted -o out.fasta

script either reformats badly formated fasta file. Within reason. Ie. Word documents will still break it.

if lenght of seq is longer than -l default 3 - writes to file.

can filter fasta by giving it a list of --not_wanted or --wanted names.

REQUIRES Biopython

"""

parser = OptionParser(usage=usage)
parser.add_option("-i", dest="in_file",
                  default=None,
                  help="current fasta you want to reformat")


parser.add_option("-l", "--length",
                  dest="length", default="5",
                  help="Output filename",
                  metavar="FILE")

parser.add_option("-o", "--out", dest="out",
                  default=None,
                  help="Output filename",
                  metavar="FILE")

(options, args) = parser.parse_args()

in_file = options.in_file
out = options.out
length = options.length

# Run as script
if __name__ == '__main__':
    reformat_as_force(in_file, length)
    reformat_as_fasta(out)
    sizes, min_contig, max_contig, \
           avg_contig, num_contig = get_fasta_stats(in_file)
    plot_multi_histogram_graph(sizes, out)
    f_out = open("gene.stats.txt", "w")
    data_out = "\t".join(["#min_contig",
                          "max_contig",
                          "avg_contig"])
    f_out.write(data_out + "\n")
    data_out = "\t".join([str(min_contig),
                          str(max_contig),
                          "%0.2f" %(avg_contig)])
    f_out.write(data_out + "\n")
    f_out.close()
    os.remove("temp_fa.fa")

