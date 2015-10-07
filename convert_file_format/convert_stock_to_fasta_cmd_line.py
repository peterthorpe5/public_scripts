from Bio import SeqIO
import os
from sys import stdin,argv
import sys
from optparse import OptionParser


script_dir = os.path.dirname(os.path.abspath(__file__))
dest_dir = os.path.join(script_dir, 'fasta_files')
try:
    os.makedirs(dest_dir)
except OSError:
    print "already exists"


def convert_files(indirectory, filename_endswith):
    for filename in os.listdir(indirectory):
        if filename.endswith(filename_endswith):
            #print filename
            #full_path_to_file = "%s/%s" (indirectory, filename)
            #full_path_to_file = os.path.join(filename)
            alignment = AlignIO.read(open(filename), "stockholm")
            outfile = "./fasta_files/%s.fasta" %(filename[:-4])

            AlignIO.write(alignment, outfile, "fasta")
    return "finished"



if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

converts

$ python convert_stock_to_fasta.py -d directory_in (default current working directory) -e file_endswith_name (or the whole filename of interest)

make sure you run this in the folder where the file is...


walk through a folder (-d) and anything that endswith -e will be converted to a fasta file

or to do one file. make -e your specific file
"""

parser = OptionParser(usage=usage)

parser.add_option("-d", dest="directory", default=default=os.getcwd(),
                  help="directory of the files you are wanting to go through")
parser.add_option("-e", "--ends", dest="endswith", default=".sto",
                  help="Output filename",
                  metavar="FILE")




(options, args) = parser.parse_args()

directory_in = options.directory
file_endswith = options.endswith

convert_files(directory_in, file_endswith)
