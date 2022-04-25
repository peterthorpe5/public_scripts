import os
from sys import stdin,argv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def seq_getter(filename1, outfile):
    "make a bed file of sequences prsent."
    f_out = open(outfile, "w")
    for seq_record in SeqIO.parse(filename1, "fasta"):
        data = "%s\t1\t%d\n" % (seq_record.id,
                                len(seq_record.seq))
        f_out.write(data)
    f_out.close()


print("""usage: python make...py fasta_file outfile

      make a bed file.
      chromosome start end""")

seq_getter(argv[1],argv[2])


