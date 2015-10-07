
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
import os
from sys import stdin,argv
import sys
from optparse import OptionParser


script_dir = os.path.dirname(os.path.abspath(__file__))
dest_dir = os.path.join(script_dir, 'Phylip_relaxed_files')
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
            alignment = AlignIO.read(open(filename), "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))
            outfile = "./Phylip_relaxed_files/%s.phy" %(filename[:-6])

            AlignIO.write(alignment, outfile, "phylip-relaxed")
    return "finished"




#this file converts an aligned fasta file and converts to a phylip  file
#alignment = AlignIO.read(open("G_B_M_seed_chopped003_more_muscle.fasta"), "fasta", alphabet=Gapped(IUPAC.protein))
#g = open("G_B_M_seed_chopped003_more_muscle.phy", "w")
#g.write (alignment.format("phylip-relaxed"))



if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

converts

$ python convert_fasta_phylip_cmd_line.py -d directory_in (default current working directory) -e file_endswith_name

make sure you run this in the folder where the file is...


walk through a folder (-d) and anything that endswith -e will be converted to a phylip file

or to do one file. make -e your specific file
"""

parser = OptionParser(usage=usage)

parser.add_option("-d", dest="directory", default=default=os.getcwd(),
                  help="directory of the files you are wanting to go through")
parser.add_option("-e", "--ends", dest="endswith", default="ed.fasta",
                  help="Output filename",
                  metavar="FILE")




(options, args) = parser.parse_args()

directory_in = options.directory
file_endswith = options.endswith

convert_files(directory_in, file_endswith)
