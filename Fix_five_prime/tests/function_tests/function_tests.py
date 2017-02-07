
#imports
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
import os
import subprocess
import tempfile
from collections import deque
from optparse import OptionParser
import datetime

def get_total_coverage(bam_file, outfile):
    ## Run samtools idxstats (this get the coverage for all transcripts:
    cmd = 'samtools idxstats "%s" > "%s"' % (bam_file, outfile)
    #data was saved in idxstats_filename
    return_code = os.system(cmd)
    if return_code:
        clean_up()
        sys_exit("Return code %i from command:\n%s" % (return_code, cmd))
    # creat a dictioanry to hold all the total expression values for the transcripts.
    overall_expression_dic = dict()
    with open(outfile, "r") as handle:
        for line in handle:
            data = line.rstrip("\n").split("\t")
            transcript = data[0]
            overall_expression_dic[transcript] = data[1:]
    #print overall_expression_dic["Mp_O_20647_c0_seq2"]
    #returns a dictionary: key[transcript], vals = ['577', '274', '0'] len, reads_mapped, last_coloumn
    return overall_expression_dic






#############################################################################

usage = """Use as
"""


parser = OptionParser(usage=usage)

parser.add_option("-t", "--transcriptome", dest="transcriptome", default=None,
                  help="the transcriptome assembly .fasta",
                  metavar="FILE")

parser.add_option("--cds", dest="cds", default=None,
                  help="the predicted cds from the transcriptome assembly .fasta",
                  metavar="FILE")

parser.add_option("--gff", dest="gff", default=None,
                  help="the predicted coordinate for the cds predictions .gff",
                  metavar="FILE")

parser.add_option("--prot", dest="prot", default=None,
                  help="the predicted amino acid cds  .pep",
                  metavar="FILE")

parser.add_option("-g", "--genome", dest="genome", default=None,
                  help="the genome sequence. Not currently used. TO DO",
                  metavar="FILE")

parser.add_option("--bam", dest="bam", default=None,
                  help="the sorted, indexed bam file as a result of the reads being"
                  "mapped back to the transcriptome  .bam",
                  metavar="FILE")
parser.add_option("--min_read_depth", dest="min_read_depth", default="10",
                  help="the min_read_depth that a transcript must have before it is"
                  "considered for statistical analysis. Default = 10")

parser.add_option("--help_full", dest="help_full", default=None,
                  help="prints out a full description of this program")

parser.add_option("-o", "--out", dest="outfile", default="results.out",
                  help="Output filename",
                  metavar="FILE")

(options, args) = parser.parse_args()


#-g
genome = options.genome
#-t
transcriptome = options.transcriptome
#--cds
cds_file = options.cds
#--gff
gff = options.gff
#--bam
bam = options.bam
#-o
outfile= options.outfile
# -min_read_depth
min_read_depth = options.min_read_depth
                
get_total_coverage(bam, "temp_out")



