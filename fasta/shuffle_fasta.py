#script to re-write badly formatted fasta file. Remove duplicates,
#or get seq or interest.

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
import os
from optparse import OptionParser
import random

    

def reformat_as_fasta(filename, outfile):
    "this function re-write a file as a fasta file"
    f = open(outfile, 'w')

    for seq_record in SeqIO.parse(filename, "fasta"):
        seq_list = []
        seq = str(seq_record.seq)
        for AA in seq:
            seq_list.append(AA)
        random.shuffle(seq_list)
        seq = ''.join(seq_list)
        seq_record.seq = seq

        seq_record.id = "shuffled_%s" % seq_record.id

        record = SeqRecord(
        Seq(seq),
        id=seq_record.id,
        description=seq_record.description,)

        SeqIO.write(record, f, "fasta") 
        
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
                  default="COVID_v_uniprot_merged.fasta",
                  help="current fasta you want to shuffle")


parser.add_option("-o", "--out", dest="out",
                  default="COVID_v_uniprot_SHUFFLED.fasta",
                  help="Output filename",
                  metavar="FILE")




(options, args) = parser.parse_args()

in_file = options.in_file
out = options.out


reformat_as_fasta(in_file, out)

