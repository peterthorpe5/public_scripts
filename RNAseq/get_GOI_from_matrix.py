# script to return gene counts for GOI
# author Pete Thorpe 2021 Oct
# imports

import os
import sys
from optparse import OptionParser
import datetime
import logging
import logging.handlers
import time


def test_line(line):
    """returns true lines. Not comments or blank line"""
    if not line.strip():
        return False  # if the last line is blank
    if line.startswith("#"):
        return False  # comment line
    return line


def parse_goi(goi, logger):
    """func take in a file with a list og gene names. 
    returns a set"""
    f_in = open(goi, "r")
    input_count = 0
    goi_set = set([])
    for line in f_in:
        input_count = input_count + 1
        if test_line(line):
            goi_set.add(line.strip())
    info = "input GOI list was %d genes" % input_count
    logger.info(info)
    info = "non-redundant GOI list was %d genes" % len(goi_set)
    logger.info(info)
    return(goi_set)


def parse_matrix(matrix, GOI_set, outfile, logger):
    """fucn take in the counts.matrix and GOI_set
    if the gene in the matrix is in the set,
    prints it to the file. """
    f_in = open(matrix, "r")
    f_out = open(outfile, "w")
    matrix_found_count = 0
    for line in f_in:
        
        if test_line(line):
            if line.startswith("\t"):
                f_out.write(line)
                continue
            data = line.split()
            name = data[0]
            if name in GOI_set:
                matrix_found_count = matrix_found_count + 1
                f_out.write(line)
    info = "GOI found in matrix was %d genes" % matrix_found_count
    logger.info(info)



if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

$ python filter....py -m matrix --goi list_of_gene_of_interest

"""

parser = OptionParser(usage=usage)

parser.add_option("-m", dest="matrix",
                  default="HLphyDis3.counts.matrix.isoform.counts.matrix",
                  help="counts.matrix")

parser.add_option("--goi", dest="goi",
                  default="GOI.txt",
                  help="file with a list of gene names",
                  metavar="FILE")
                  
parser.add_option("-o", dest="outfile",
                  default="GOI.counts.matrix",
                  help="out file.matrix")



(options, args) = parser.parse_args()

matrix = options.matrix
goi = options.goi
outfile = options.outfile

logfile = "GOI_counts.log" 
# Run as script
if __name__ == '__main__':
    start_time = time.time()
    # Set up logging
    logger = logging.getLogger('get_GOI.py: %s'
                               % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)
    try:
        logstream = open(logfile, 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        logger.error("Could not open %s for logging" %
                     logfile)
        sys.exit(1)
    # get gene of interests
    GOI_set = parse_goi(goi, logger)
    # parse the matrix 
    parse_matrix(matrix, GOI_set, outfile, logger)
            
