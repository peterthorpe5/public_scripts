#!/usr/bin/env python3
# title: reduce the info in the output of SNPeff
# the output is useful of course, but stops it being readble. 
# Author: Peter Thorpe. 


# imports
import sys
import os
from optparse import OptionParser
import datetime
import logging
import logging.handlers
import time
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def index_genome_file(genome, logger):
    """index the genome file
    return a seq_record object that can be accessed
    in a dictionary like manner"""
    logger.info("indexing genome: %s", genome)
    genome_database = SeqIO.index(genome, "fasta")
    logger.info("indexing genome finished")
    return genome_database


def parse_gff(gff, genome, outfile, logger):
    """take in gff file. """
    f_in = open(gff, "r")
    f_out = open(outfile, "w")
    genome_database = index_genome_file(genome, logger)
    last_scaff = ""
    scaff_count = 0
    print_count = 0
    scaff_list = set([])
    for line in f_in:
        if line.startswith("#"):
            f_out.write(line)
            continue
        if not line.strip():
            continue # end of file
        assert len(line.split("\t")) == 9 , "GFF fields wrong length should be 9"
        scaff, source, feature, start, stop, score, \
              direction, frame, gene_info = line.split("\t")
        scaff_list.add(scaff)
        if last_scaff == "":
            last_scaff = scaff
            scaff_count = scaff_count + 1
        if scaff == last_scaff:
            f_out.write(line)
        if scaff != last_scaff:
            last_scaff = scaff
            print_count = print_count + 1
            scaff_count = scaff_count + 1
            f_out.write(line)
    f_out.write("##FASTA\n")
    logger.info("writting this many record %d", len(scaff_list))
    for entry in sorted(scaff_list):
        seq_record = genome_database[entry.rstrip()]
        seq_record.description= ""
        # logger.info("writting record %s", entry)
        SeqIO.write(seq_record, f_out, "fasta")
    logger.info("number scaff = %d, we printed %d", scaff_count, print_count)
    f_out.close()
    f_in.close()
            
            
        
        

usage = """ python make ....py -h

python reduce_gff.py --gff snpef.gff -o outreduced.gff


# the output is useful of course, but stops it being readble. 
# Author: Peter Thorpe. 

"""



parser = OptionParser(usage=usage)


parser.add_option("--gff", dest="gff",
                  default="temp.gff",
                  help="gff",
                  metavar="FILE")

parser.add_option("--logger", dest="logger",
                  default=None,
                  help="Output logger filename. Default: " +
                  "outfile_std.log",
                  metavar="FILE")

parser.add_option("--genome", dest="genome",
                  default=None,
                  help="genome")

parser.add_option("-o", "--out",
                  dest="outfile",
                  default="combined_gff_fa.gff",
                  help="combined_gff_fa.gff",
                  metavar="FILE")


(options, args) = parser.parse_args()



#--gff
gff = options.gff
# genome
genome = options.genome
# outfile
outfile = options.outfile



#######################################################################
# Run as script
# Run as script
if __name__ == '__main__':
    # Set up logging
    if not options.logger:
        options.logger = "reduceSnpEff_gff.log"
    logger = logging.getLogger('TranStart.py: %s' % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)
    try:
        logstream = open(options.logger, 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        outstr = "Could not open %s for logging" % options.logger
        logger.error(outstr)
        sys.exit(1)
    # Report input arguments
    logger.info(sys.version_info)
    logger.info("Command-line: %s", ' '.join(sys.argv))
    logger.info("Starting testing: %s", time.asctime())
    if not os.path.isfile(gff):
        print("gff not here motherfucker!")
    parse_gff(gff, genome, outfile, logger)
