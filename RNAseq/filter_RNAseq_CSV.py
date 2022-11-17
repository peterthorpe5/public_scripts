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
import pandas as pd



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
    df = pd.read_csv(infile)
    df = df.reset_index()
    f_out = open(outfile, "w")
    gene_names = set([])
    for index, row in df.iterrows():
        # incositant naming in files, so use pandas
        name = row["row"]
        logfc = row["log2FoldChange"]
        FDR = row["padj"]
        out_data = "%s\t%s\t%s\n" % (name, logfc, FDR)
        if name in name_set:
            # print("duplicate", name)
            continue
        # add the name to the set        
        name_set.add(name)
        if logfc >= float(LOGFC_threshold):
            #print(line.rstrip())
            if FDR <= float(FDR_threshold):
                count = count + 1
                gene_names.add(name)
                f_out.write(out_data)
        # convert the negative to positive for easy testing
        if (logfc * -1.0) >= float(LOGFC_threshold):
            if FDR <= float(FDR_threshold):
                count = count + 1
                gene_names.add(name)
                out_data = "%s\t%s\t%s\n" % (name, logfc, FDR)
                f_out.write(out_data)

    f_out.close()
    gene_names = "\t".join(gene_names)
    data =("%s\twe found %d DE transcripts\ttranscript:\t%s" % (infile,
                                                                count,
                                                                gene_names))
    logger.info(data)



if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

$ python filter....py --fdr FDR_threshold --log logFC threshold

"""

parser = OptionParser(usage=usage)

parser.add_option("--fdr", dest="fdr",
                  default=0.01,
                  help="FDR threshold: default =< 0.01")


parser.add_option("-l", "--logfc", dest="out",
                  default=1.5,
                  help="logFC threshold: default => 1.5",
                  metavar="FILE")



(options, args) = parser.parse_args()

FDR = options.fdr
LOGFC = options.out

FDR = float(FDR)
LOGFC = float(LOGFC)

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
        if not filename.endswith("_allresults.csv") : continue
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
        str_LOGFC = "%.2f" % LOGFC
        outfile = "%s_LOGFC_%s_FDR_%s" % (infile.split("_allresults.csv")[0],
                                         str_LOGFC,
                                         str(FDR))
        outfile = os.path.join(out_folder, outfile)
        parse_file(infile, outfile, LOGFC, FDR, logger)
            
