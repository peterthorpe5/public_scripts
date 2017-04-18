#script to re-write badly formatted fasta file. Remove duplicates,
#or get seq or interest.

import sys
import os
from optparse import OptionParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


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


def reformat_as_fasta(outfile):
    "this function re-write a file as a fasta file"
    f= open(outfile, 'w')
    for seq_record in SeqIO.parse("temp_fa.fa", "fasta"):
        SeqIO.write(seq_record, f, "fasta")
    f.close()
    return True


def reformat_as_fasta(filename, length):
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
    return True



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
                  dest="length", default="20",
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
    reformat_as_fasta(in_file, length)
    reformat_as_fasta("temp_fa.fa", out)
    min_contig, max_contig, avg_contig, num_contig = get_fasta_stats(out)
    f_out = open("gene.stats.txt", "w")
    data_out = "\t".join(["#min_contig",
                          "max_contig",
                          "avg_contig"])
    f_out.write(data_out + "\n")
    data_out = "\t".join([str(min_contig),
                          str(max_contig),
                          str(avg_contig)])
    f_out.write(data_out + "\n")
    f_out.close()
    os.remove("temp_fa.fa")

