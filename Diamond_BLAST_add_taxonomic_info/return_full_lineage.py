#!/usr/bin/env python

# author: Peter Thorpe September 2020. U of St A, Dundee, UK.

# imports
from __future__ import print_function
import time
import os
import sys
from optparse import OptionParser  # TODO: update to argparser
import datetime
import logging
import logging.handlers

# get the files:
# wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
# tar -zxvf new_taxdump.tar.gz


def test_line(line):
    """returns true lines. Not comments or blank line"""
    if not line.strip():
        return False  # if the last line is blank
    if line.startswith("#"):
        return False  # comment line
    return line



def tax_id_to_lineage(rankedlineage):
    """function to return a dict of tax_id to lineage
    # wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
    # tar -zxvf new_taxdump.tar.gz
    file is a dark side format: 12199   |       Cowpea aphid-borne mosaic virus |     

    """
    tax_id_to_lineage = dict()
    with open(rankedlineage, "r") as handle:
        for line in handle:
            if not test_line(line):
                continue
            fields = [x.strip() for x in line.rstrip("\t|\n").split("\t|\t")]
            taxid = int(fields[0])
            full_lineage = fields[1:]
            # print(fields)
            reverse_order = "| ".join(full_lineage[::-1])
           
            tax_id_to_lineage[str(taxid)] = reverse_order
    return tax_id_to_lineage


def tax_names(tax_id_file, tax_id_to_lineage,
              outfile, logger):
    """read in the tab file. 1 tax id per line.
    Queries the dictionary for the name and report the
    linegae in reverse order back. 
    """
    print(dict(list(tax_id_to_lineage.items())[0:2]))
    f_out = open(outfile, "w")
    with open(tax_id_file) as file:
        for line in file:
            if not test_line(line):
                continue
            tax_id = line.rstrip()
            if tax_id in tax_id_to_lineage:
                lineage = tax_id_to_lineage[tax_id]
                f_out.write(tax_id + "\t" + lineage + "\n")
            else:
                logger.warning("%s has no lineage info", tax_id)
                f_out.write(tax_id + "\t" + "NA" + "\n")
    f_out.close()
            
                                     

#############################################################################


if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)

usage = """Use as follows:
$ python return_full_lineage.py -i list_of_tax_ids.txt -p path_to_tax_files -o outfile


    or

$ python return_full_lineage.py -i list_of_tax_ids.txt -r rankedlineage.dmp -o outfile


# tax_id file.

1 tax id per line. You can use linux cut to get the coloumn of interest from your
tax_id_annotated blast out put.

cat -f(what_ever_column_it_is) > tax_id_of_interest.

now you have your blast file and a corresponding lineage for the hits.
If you are good with awk you can merge them. 

"""

parser = OptionParser(usage=usage)

parser.add_option("-i", "--in",
                  dest="tax_id_file",
                  default="wanted",
                  help="file containing 1 tax id per line " +
                  "cut the blast output to get this, if required",
                  metavar="FILE")

parser.add_option("-p", "--path",
                  dest="path",
                  default=os.getcwd(),
                  help="Directory containing taxonomy/database files " +
                  "(set by -t, -c, -n, -d). Default = current working " +
                  "directory. This is not used with the main input and output "
                  "filenames. " +
                  "Dir = where you put all the " +
                  "downloaded NCBI files................... IF YOU GIVE " +
                  "THE PATH YOU DONT NEED TO SET -t, -c, -n, -d))")

parser.add_option("-r", "--rankedlineage",
                  dest="rankedlineage",
                  default="rankedlineage.dmp",
                  help="NCBI provided file prot.accession2taxid " +
                  "after unzipping (from FTP site, "
                  " after unzipping). " +
                  "These file required file options can be left blank if -p " +
                  "is specified with a path to where all these can be found. " +
                  "If -p /PATH/ is specified python will look in the " +
                  "folder by default.",
                  metavar="FILE")

parser.add_option("-o", "--out",
                  dest="outfile",
                  default="full_lineage.txt",
                  help="Output filename - " +
                  "",
                  metavar="FILE")

(options, args) = parser.parse_args()


def apply_path(folder, filename):
    """If filename is not absolute, assumed relative to given folder.

    Here filename is a relative path (does not start with slash):

    >>> apply_path("/mnt/shared/taxonomy", "names.dmp")
    '/mnt/shared/taxonomy/names.dmp'

    Here filename is already an absolute path, so no changes:

    >>> apply_path("/mnt/shared/taxonomy", "/tmp/ncbi/taxonomy/names.dmp")
    '/tmp/ncbi/taxonomy/names.dmp'

    """
    if os.path.isabs(filename):
        return filename
    else:
        return os.path.join(folder, filename)

# -i
tax_id_file = options.tax_id_file
#-r
rankedlineage = apply_path(options.path, options.rankedlineage)
# -p
path_files = options.path
#-o
outfile = options.outfile


# Run as script
if __name__ == '__main__':
    # Set up logging
    log_out = outfile + ".log"
    logger = logging.getLogger('ranked_lineage.py: %s' % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)
    try:
        logstream = open(log_out, 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        print("Could not open %s for logging" % log_out)
        sys.exit(1)
    # Report input arguments
    logger.info(sys.version_info)
    logger.info("Command-line: %s", ' '.join(sys.argv))
    logger.info("Starting testing: %s", time.asctime())
    # call the main function
    filename_list = [rankedlineage,
                     tax_id_file]
    for needed_file in filename_list:
        if not os.path.isfile(needed_file):
            logger.info("sorry cannot find you %s file", needed_file)
            logger.info("please check this command again, " +
                        "with the full path if required")
            os._exit(0)
    tax_id_to_lineage = tax_id_to_lineage(rankedlineage)
    tax_names(tax_id_file, tax_id_to_lineage, outfile, logger)

