
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
#os imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser



def reformat_as_fasta(filename,length,outfile):
    "this function re-write a file as a fasta file"
    f= open(outfile, 'w')

        
    #print wanted_data
    for seq_record in SeqIO.parse(filename, "fasta"):
        seq_record.id = str(seq_record.id).split("_")[0]
        seq_record.id= seq_record.id.replace("sca", "RpS")
        seq_record.description = ""
        SeqIO.write(seq_record, f, "fasta")                    
    
    f.close()
    return True



if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

converts

$ python rewrite_as_fasta.py -i in.fasta -l min_length_of_seq (default(3)) --not_wanted --wanted -o out.fasta

script either reformats badly formated fasta file. Within reason. Ie. Word documents will still break it.

if lenght of seq is longer than -l default 3 - writes to file.

can filter fasta by giving it a list of --not_wanted or --wanted names.

"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="in_file", default=None,
                  help="current fasta you want to reformat")


parser.add_option("-l", "--lenth", dest="length", default="3",
                  help="Output filename",
                  metavar="FILE")
parser.add_option("-o", "--out", dest="out", default=None,
                  help="Output filename",
                  metavar="FILE")




(options, args) = parser.parse_args()

in_file = options.in_file

out = options.out
length = options.length


reformat_as_fasta(in_file,length,out)
print 'done'

