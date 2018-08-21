# Example parsing/processing of Phytophthora FASTA read files
import math
import os
import argparse
# NOTE TO SELF: NORMALISE AGAINST STARING NUMBER OF READS. 
from collections import defaultdict
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
from Bio import SeqIO
import sys


VERSION = "summerise results: v0.01"
if "--version" in sys.argv:
    print(VERSION)
    sys.exit(1)

usage = """

"""
if "--help" or "-h" in sys.argv:
    print(usage)

# Generate counts of unique sequences for each input file
count_dict = defaultdict(int)
counts_by_file = {}
#fnames = ("file1.fasta", "file2.fasta")

def get_args():
    parser = argparse.ArgumentParser(description=" count_seq_abundance_from_related_samples.py",
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("count_seq_abundance_from_related_samples")[0]
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument("-i", "--in", dest='in_files',
                          action="store",
                          default=None,
                          type=str,
                          help="in_files, commor seperated list")
    optional.add_argument("-o", "--out", dest='out',
                          action="store",
                          default="NO",
                          type=str,
                          help="outfilename")
    optional.add_argument("-s", "--seq", dest='seq',
                          action="store",
                          default="NO",
                          type=str,
                          help="print unquie sequences")
    optional.add_argument('--version',
                          action='version',
                          version="%s: abundance.py " + VERSION)
    args = parser.parse_args()
    return args, file_directory

out_names = ""
args, FILE_DIRECTORY = get_args()
if "," in args.in_files:
    files = args.in_files.split(",")
    for infile in files:
        name = os.path.split(infile)[-1]
        name = name.split("_L001.assembled.fastq")[0]
        # print(name)
        out_names = out_names + name + "_VS_"
    out_names = out_names.rstrip("_VS_")
else:
    files = [args.in_files]
    out_names = os.path.split(files)[-1].split("_L001.assembled.fastq")[0]
for fname in files:
    #print("files", files)
    #print("fname =", fname)
    infile = os.path.split(fname)[-1]
    #print("infile = ", infile)
    fstem = os.path.splitext(infile)[0]
    counts_by_file[fstem] = defaultdict(int)
    for record in SeqIO.parse(fname, "fasta"):
        count_dict[str(record.seq)] += 1
        counts_by_file[fstem][str(record.seq)] += 1

# Summary info:
print("Total unique read sequences: {}".format(len(count_dict)))
for key, val in counts_by_file.items():
    print("Total unique read sequences in {}: {}".format(key, len(val)))

# Show sorted list of counts of reads from each file
for name, cdict in counts_by_file.items():
    print("Working with file: {}".format(name))
    data = sorted(cdict.items(), key=lambda x: x[1], reverse=True)
    #print("\n\t"+ "\n\t".join(["{}: {}".format(key, val) for key, val in data[:10]])+ "\n")

# How many sequences in common
common_seqs = set(count_dict.keys())
print("Total unique sequences: {}".format(len(common_seqs)))
for key, val in counts_by_file.items():
    print("Working with sequences from {}: {}".format(key, len(val)))
    print("Finding intersection...")
    common_seqs.intersection_update(val.keys())
    print("New count of common sequences: {}".format(len(common_seqs)))

# if option to shwo seq:
if args.seq.upper != "NO":
    count = 0
    for entry in common_seqs:
        count += 1
        print ">common_seq_%d\n%s" % (count, entry)

# Report on common sequences
common_dict = {key: value for key, value in count_dict.items() if key in common_seqs}
data = sorted(common_dict.items(), key=lambda x: x[1], reverse=True)
#print("\n\t" + "\n\t".join(["{}: {}".format(key, val) for key, val in data[:100]]) + "\n")

# Plot common counts
plt.style.use('seaborn-darkgrid')

plotdata = enumerate([item[1] for item in data])
x, y = list(zip(*plotdata))
fig = plt.plot(x, [math.log(val) for val in y], marker='o', markersize=2)
if args.out.upper() != "NO":
    outname = args.out
else:
    outname = out_names
plt.grid(True)
plt.title("Data: %s" % (outname),
              fontsize=6)
plt.xlabel("Number of of times a sequences occurs")
plt.ylabel("LOG Number of unique sequences")
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.xticks(fontsize=6, rotation=90)
plt.show()
# incase the name is super long, set withthe -o option 
plt.savefig(outname + "_plot" + ".png")
