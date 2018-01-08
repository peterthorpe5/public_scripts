
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
from sys import stdin,argv
import sys
from optparse import OptionParser


def reformat_as_fasta(input_handle, PREFIX, outfile):
    "this function re-write a file as a fasta file"
    output_handle = open(outfile, 'w')
    count = 0
    for seq_record in SeqIO.parse(input_handle, "genbank") :
        print "Dealing with GenBank record %s" % seq_record.id
        for seq_feature in seq_record.features :
            if seq_feature.type=="CDS" :
                count += 1
                assert len(seq_feature.qualifiers['translation'])==1
                output_handle.write(">%s%04d\t%s from %s\n%s\n" % (
                       PREFIX, count, seq_feature.qualifiers['product'],
                       seq_record.name,
                       seq_feature.qualifiers['translation'][0]))

    output_handle.close()

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

converts

$ python rewrite_as_fasta.py -i in.gbk -p prefix_for_gene_names -o out.fasta

script to extract cds from gbk files.
Output the CDS as
>prefix00001
AHGATGHD etc ...

"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="in_file", default=None,
                  help="current fasta you want to reformat")
parser.add_option("-p", "--prefix", dest="prefix", default="gene",
                  help="prefix_for_genes to be called",)
parser.add_option("-o", "--out", dest="out", default="AA.fasta,
                  help="Output filename",
                  metavar="FILE")


(options, args) = parser.parse_args()

in_file = options.in_file
prefix = options.prefix
out = options.out

if __name__ == '__main__':
    reformat_as_fasta(in_file, prefix, out)

