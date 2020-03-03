
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import generic_dna
#os imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser
from Bio.Data import CodonTable
from collections import defaultdict



def get_codon(seq, codon_count):
    """funct to get the codon on batches of 3s"""

    for i in range(3, len(seq) + 1, 3):
        codon = seq[i-3:i]
        if "N" in codon:
            return codon_count
        codon_count[codon] += 1
    return codon_count
    


def write_out(codon_count, outfile):
    """func to write out the codon dict"""
    f_out = open(outfile, "w")
    title = "#codon\ttranslation\tfrequency\n"
    f_out.write(title)
    for codon, counts in codon_count.items():
        seq_record = Seq(codon, generic_dna)
        protien = seq_record.translate()
        outfmt = "%s\t%s\t%d\n" % (codon, protien,
                                   counts)
        f_out.write(outfmt)
    f_out.close()


def reformat_as_fasta(filename, outfile):
    "this function re-write a file as a fasta file"
    f = open(outfile, 'w')
    codon_count = defaultdict(int)

    for seq_record in SeqIO.parse(filename, "fasta"):
        seq = str(seq_record.seq.upper())
        codon_count = get_codon(seq, codon_count)
    write_out(codon_count, outfile)
    f.close()



if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

trims 20 bp off seq

$ python rewrite_as_fasta.py -i in.fasta --o out.fasta


"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="in_file",
                  default="H_sch_gene_calls_v1.codingseq",
                  help="current fasta you want to reformat")

parser.add_option("-o", "--out", dest="out",
                  default="codon_freq_table.txt",
                  help="Output filename",
                  metavar="FILE")


(options, args) = parser.parse_args()

in_file = options.in_file
out = options.out

reformat_as_fasta(in_file, out)


