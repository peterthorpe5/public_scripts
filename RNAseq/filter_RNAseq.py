# script to parse and filter DE results for sig FC and FDR
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


def parse_file(infile, outfile, LOGFC_threshold, FDR_threshold, logger):
    """func to parse the de results and writes the sig results"""
    count = 0
    name_set = set([])
    f = open(infile, "r")
    f_out = open(outfile, "w")
    for line in f:
        if test_line(line):
            if "Row.names" in line:
                f_out.write(line.lstrip())
                continue
            # as my data currently is, row are just line numbers
            line = "\t".join(line.split("\t")[1:])
            data = line.split()
            name = data[0]
            if name in name_set:
                print("duplicate", name)
                continue
            
            name_set.add(name)
            logfc = data[3]
            FDR = data[7]
            logfc =  float(logfc)
            FDR = float(FDR)
            if logfc >= float(LOGFC_threshold):
                #print(line.rstrip())
                if FDR <= float(FDR_threshold):
                    count = count + 1
                    f_out.write(line)

    f.close()
    f_out.close()
    data =("%s\twe found %d DE transcripts" % (infile, count))
    logger.info(data)



if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

$ python filter....py --fdr FDR_threshold --log logFC threshold

"""

parser = OptionParser(usage=usage)

parser.add_option("--fdr", dest="fdr",
                  default=0.05,
                  help="FDR threshold: default =< 0.05")


parser.add_option("-l", "--logfc", dest="out",
                  default=1.0,
                  help="logFC threshold: default => 1.0",
                  metavar="FILE")



(options, args) = parser.parse_args()

FDR = options.fdr
LOGFC = options.out

logfile = "filter_RNAseqogFC_%.2f_FDR_%.2f.log" % (LOGFC, FDR)
# Run as script
if __name__ == '__main__':
    start_time = time.time()
    # Set up logging
    logger = logging.getLogger('filter_RNAseq.py: %s'
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

    # iterate through all the .DE_results files.
    for filename in os.listdir(".") :
        if not filename.endswith(".DE_results") : continue
        FDR = float(FDR)
        LOGFC = float(LOGFC)
        # folder_for_all_the_out_files
        out_folder = "filtered_by_logFC_%.2f_FDR_%.2f" % (LOGFC, FDR)
        dest_dir = os.path.join(os.getcwd(),
                                out_folder)
        try:
            os.makedirs(dest_dir)
        except OSError:
            pass
        infile = filename
        outfile = "%s_LOGFC_%s_FDR_%s" % (infile.split("_results")[0],
                                         str(LOGFC),
                                         str(FDR))
        outfile = os.path.join(out_folder, outfile)
        parse_file(infile, outfile, LOGFC, FDR, logger)
            
