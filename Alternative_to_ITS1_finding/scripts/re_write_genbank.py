# script to reformat fasta to nexus format
from Bio import SeqIO
#os imports
import os


def convert_file(in_file, out_file):
    sequences = SeqIO.parse(in_file, "genbank")
    g = open(out_file, "w")
    SeqIO.write(sequences, out_file, "genbank")



for filename in os.listdir("."):
    if not filename.endswith(".gb") : continue
    convert_file(filename, filename.split(".gb")[0]+".gbk")
