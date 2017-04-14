#script to re-write badly formatted fasta file. Remove duplicates,
#or get seq or interest.

import sys
import os
from optparse import OptionParser
    

def reformat_as_fasta(filename, length, outfile):
    "this function re-write a file as a fasta file"
    f = open(outfile, 'w')
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
    reformat_as_fasta(in_file, length, out)

