#script to re-write badly formatted fasta file. Remove duplicates,
#or get seq or interest.

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
import os
from optparse import OptionParser

    

def reformat_as_fasta(filename, in2, outfile):
    "this function re-write a file as a fasta file"
    f = open(outfile, 'w')
  
    #print wanted_data
    record = []
    seq_set = set([])
    id_set = set([])
    duplicate_count = 0
    file1_count = 0
    file2_count = 0

    for seq_record in SeqIO.parse(filename, "fasta"):
        file1_count = file1_count + 1
        if str(seq_record.seq) not in seq_set:
            seq_set.add(str(seq_record.seq))
            record.append(seq_record)
            id_set.add(seq_record.id)
        else:
            print("Duplicate sequence found", seq_record)
            duplicate_count += 1

    for seq_record in SeqIO.parse(in2, "fasta"):
        file2_count = file2_count + 1

        if str(seq_record.seq) not in seq_set:
            seq_set.add(str(seq_record.seq))
            record.append(seq_record)
            id_set.add(seq_record.id)
        else:
            print("Duplicate sequence found", seq_record)
            duplicate_count += 1

    print("going to write %d record" % len(record))
    print("found %d duplicate sequneces " % duplicate_count)
    print("file1 had %d seq, file2 had %d seq: total = %d " % (file1_count,
                                                               file2_count,
                                                          file2_count + file1_count))

    print("going to write %d seqs " % (len(record)))
    out_set = set([])
    for seq_record in SeqIO.parse(filename, "fasta"):
        if seq_record.id in id_set:
            if seq_record.id not in out_set:
                out_set.add(seq_record.id)
                SeqIO.write(seq_record, f, "fasta")
    print("first file done")
    for seq_record in SeqIO.parse(in2, "fasta"):
        if seq_record.id in id_set:
            if seq_record.id not in out_set:
                out_set.add(seq_record.id)
                SeqIO.write(seq_record, f, "fasta") 
        
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
                  default="NC_045512.txt",
                  help="current fasta you want to merge")


parser.add_option("-p", dest="in2",
                  default="uniprot-filtered-reviewed_yes+AND+organism__Homo+sapiens+(Human)+[96--.fasta",
                  help="second file to merge")

parser.add_option("-o", "--out", dest="out",
                  default="COVID_v_uniprot_merged.fasta",
                  help="Output filename",
                  metavar="FILE")




(options, args) = parser.parse_args()

in_file = options.in_file
in2 = options.in2
out = options.out


reformat_as_fasta(in_file, in2, out)

