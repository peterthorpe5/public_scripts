#!/usr/bin/env python
#
# metapy_sanger_read.py
#
# Code for script to identify OTUs from metabarcoding reads.
# runs multiple clustering programs
# THIS IS FOR SANGER READS
#
# (c) The James Hutton Institute 2016-2017
# Author: Peter Thorpe

import sys
import shutil
import errno
import logging
import logging.handlers
import multiprocessing
import os
import time
import traceback
import subprocess
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from pycits.tools import convert_fq_to_fa, NotExecutableError, trim_seq

from pycits.metapy_tools import decompress, compress,\
     last_exception, metapy_trim_seq, covert_chop_read, make_folder,\
     test_reads_exist_and_suffix, database_checker,\
     get_sizes, db_len_assembled_len_ok, stats_on_list_of_sizes,\
     plot_seq_len_histograms

from pycits.metapy_tools_sanger import convert_ab1_to_fq, \
     sanger_extra_qc_trim, plot_trace

from pycits import tools, \
     clean_up, muscle, trimmomatic

if sys.version_info[:2] != (3, 5):
    # e.g. sys.version_info(major=3, minor=5, micro=2,
    # releaselevel='final', serial=0)
    # break the program
    if sys.version_info[:2] != (2, 7):
        print ("currently using:", sys.version_info,
               "  version of python")
        raise ImportError("Python 3.5 or 2.7 is required for " +
                          "metapy_sanger_read.py")
        print ("did you activate the virtual environment?")
        sys.exit(1)

VERSION = "Pycits classify OTU using Sanger ab1 files: v1.0.0"
if "--version" in sys.argv:
    print(VERSION)
    sys.exit(1)

#########################################################################


def get_args():
    parser = argparse.ArgumentParser(description="Pipeline: cluster " +
                                     "data for metabarcoding ",
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("metapy")[0]
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--thread", dest='threads',
                          action="store", default="8",
                          type=str,
                          help="number of threads")

    optional.add_argument("-a", "--ab1", dest='ab1',
                          action="store",
                          default=os.path.join(file_directory,
                                               "sanger_read",
                                               "SQ16_Example.ab1"),
                          type=str,
                          help="sanger ab1 file")

    optional.add_argument("-d", "--OTU_DB", dest='OTU_DB',
                          action="store",
                          default="/mnt/scratch/local/blast/ncbi/nt",
                          type=str,
                          help="database of seq of to compare against")

    optional.add_argument("-e", "--evalue", dest='evalue',
                          action="store",
                          default=1e-05,
                          type=float,
                          help="evalue to filter results with")

    optional.add_argument("-m", "--mismatches", dest='mismatches',
                          action="store",
                          default="3",
                          type=str,
                          help="number of mismatches to filter results with")

    optional.add_argument("--left_trim", dest='left_trim',
                          action="store",
                          default=33,
                          type=int,
                          help="left_trim for primers or conserved " +
                          "regions. Default 25 ")

    optional.add_argument("--right_trim", dest='right_trim',
                          action="store",
                          default=150,
                          type=int,
                          help="right_trim for primers or conserved " +
                          "regions.  Default 0")

    optional.add_argument("--adaptors", dest='adaptors',
                          action="store",
                          default=os.path.join(file_directory,
                                               "adapters",
                                               "TruSeq3-PE.fa"),
                          type=str,
                          help="adaptors for trimming. Can supply custom " +
                          " file if desired")

    optional.add_argument("--phred", dest='phred',
                          action="store",
                          default="phred33",
                          type=str,
                          help="phred33 is default. " +
                          "Dont change unless sure")

    optional.add_argument("--verbose", dest="verbose",
                          action="store_true",
                          default=False,
                          help="Report verbose output")

    optional.add_argument("--align", dest="align",
                          action="store_true",
                          default=True,
                          help="to align clusters in the output " +
                          "you must have muscle in your PATH as muscle")


    optional.add_argument("--muscle", dest="muscle",
                          action="store", default="muscle",
                          help="Path to MUSCLE... If version already" +
                          "in PATH then leave blank")


    optional.add_argument("--trimmomatic", dest="trimmomatic",
                          action="store", default="trimmomatic",
                          help="Path to trimmomatic... If version already" +
                          "in PATH then leave blank")

    optional.add_argument("--logfile",
                          dest="logfile",
                          action="store",
                          default="pipeline.log",
                          type=str,
                          help="Logfile name")

    optional.add_argument("--cleanup",
                          dest="cleanup",
                          action="store",
                          default="yes",
                          help="deletes most files the program creates ")

    optional.add_argument("-h", "--help",
                          action="help",
                          default=argparse.SUPPRESS,
                          help="Displays this help message"
                          " type --version for version")
    optional.add_argument('--version',
                          action='version',
                          version="%s: metapy.py " + VERSION)
    args = parser.parse_args()
    return args, file_directory

###################################################################
# Global variables
WARNINGS = ""
args, FILE_DIRECTORY = get_args()
# setting up some test variables
THREADS = args.threads
READ_PREFIX = os.path.split(args.ab1)[-1].split(".ab1")[0]
PHREDSCORE = args.phred
PREFIX = READ_PREFIX
ADAPTERS = args.adaptors
OUTDIR_TRIM = "trimmed_reads"
OUTFILES = [os.path.join(OUTDIR_TRIM, PREFIX + suffix) for suffix in
            ("_paired_R1.fq.gz", "_unpaired_R1.fq.gz",
             "_paired_R2.fq.gz", "_unpaired_R2.fq.gz")]

OTU_DATABASE = args.OTU_DB
WORKING_DIR = os.getcwd()

# left_trim 53 - this should be default??
LEFT_TRIM = args.left_trim
RIGHT_TRIM = args.right_trim
RESULT_FOLDER = PREFIX + "_RESULTS"
make_folder(RESULT_FOLDER, WORKING_DIR)
RESULTS = []
#####################################################################

# TODO: There is *a lot* of repeated code here.
#       The many try-except structures
#       should be making you think that you need a function

def check_tools_exist(WARNINGS):
    """function to check to see what tools are in the PATH,
    decide what we can use
    Returns a list of programs that were exectable and a warning string.
    The warnings are tools that were not executable"""
    tools_list = []
    Warning_out = WARNINGS + "Tool executable warning: "

    try:
        trimmomatic.Trimmomatic(args.trimmomatic)
        tools_list.append("trimmomatic")
    except ValueError:
        Warning_out = Warning_out + "trimmomatic not in path\n"

    try:
        muscle.Muscle(args.muscle)
        tools_list.append("muscle")
    except ValueError:
        Warning_out = Warning_out + "muscle not in path\n"
    return tools_list, Warning_out


def make_parse_cmd(outfile, mismatches):
    """function to return the parsing command"""
    cmd_parse = " ".join(["python",
                          os.path.join(FILE_DIRECTORY,
                                       "bin",
                                       "BLAST_xml_parser.py"),
                          "-i",
                          xml_out,
                          "-o",
                          outfile
                          ,
                          "-e",
                          str(args.evalue),
                          "-m",
                          str(mismatches)])
    return cmd_parse


def run_parse_cmd(cmd, logger):
    """func to run the parser command as this may be called a few times"""
    logger.info("%s make cmd_parse command", cmd)
    #  removed check-True
    pipe = subprocess.run(cmd_parse, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    # logger.info("XML parser stdout: %s", pipe.stdout)
    # logger.info("XML parser stderr: %s", pipe.stderr)


#######################################################################
# Run as script
if __name__ == '__main__':
    # Set up logging
    logger = logging.getLogger('metapy.py: %s' % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)
    if args.logfile == "pipeline.log":
        args.logfile = PREFIX + "_pipeline.log"
    try:
        logstream = open(args.logfile, 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        outstr = "Could not open %s for logging" % args.logfile
        logger.error(outstr)
        sys.exit(1)
    # Report input arguments
    logger.info(sys.version_info)
    logger.info("Command-line: %s", ' '.join(sys.argv))
    logger.info("Starting testing: %s", time.asctime())
    logger.info("using database: %s", OTU_DATABASE)
    # Get a list of tools in path!
    logger.info("checking which programs are in PATH")
    tools_list, Warning_out = check_tools_exist(WARNINGS)
    logger.info(Warning_out)

    # convert the ab1 file to fastq then fasta
    fq_out = os.path.join(RESULT_FOLDER,  READ_PREFIX + ".fastq")
    fa_out_all = os.path.join(RESULT_FOLDER,
                              READ_PREFIX + "all.fasta")
    fa_out = os.path.join(RESULT_FOLDER,  READ_PREFIX + ".fasta")
    fa_out_QC = os.path.join(RESULT_FOLDER,  READ_PREFIX + "qc.fasta")
    convert_ab1_to_fq(args.ab1, fq_out)
    abi_file = os.path.split(args.ab1)[-1]
    image_file_out = os.path.join(RESULT_FOLDER, abi_file + ".png")
    plot_trace(args.ab1, image_file_out)

    # tests quality score are !!!! which is Q0
    fq_quality_values = "Yes"
    with open(fq_out) as f_handle:
        line_count = 0
        for line in f_handle:
           line_count = line_count + 1
           if line_count == 4:
               if line.startswith("!!!!!!!!!!!!!!!"):
                   logger.warning("Sanger files from JHI have no " +
                                  "quality values. A chopping " +
                                  "approach will be used instead. " +
                                  "This is NOT ideal!!!")
                   fq_quality_values = "No"
                   logger.info("converting fq to fa")
                   convert_fq_to_fa(fq_out, fa_out_all)
                   # metapy_trim_seq(infname, outfname,
                   #                 lclip=53, rclip=0, minlen=100)
                   metapy_trim_seq(fa_out_all, fa_out, args.left_trim,
                                   args.right_trim)
                   logger.info("QC trimming")
                   sanger_extra_qc_trim(fa_out_all, fa_out_QC, 3, logger)

    ####################################################################
    # trimmomatic trm reads
    if fq_quality_values == "Yes":
        if "trimmomatic" in tools_list:
            trim_out = os.path.join(RESULT_FOLDER, "trimmed.fq")
            trim = " ".join(["trimmomatic",
                             "SE -phred33",
                             fq_out,
                             trim_out,
                             "ILLUMINACLIP:TruSeq3-SE:2:30:10",
                             "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:5",
                             "MINLEN:36"])
            logger.info("trimmomatic command  = %s " % trim)
            pipe = subprocess.run(trim, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True)

            convert_fq_to_fa(trim_out, fa_out_QC)
            os.remove(trim_out)
            logger.warning("deleting: %s", trim_out)
    ####################################################################
    # blast, make blast db
    db_name = os.path.split(OTU_DATABASE)[-1]
    if not os.path.isfile(OTU_DATABASE + ".nhr"):
        if not OTU_DATABASE == "/mnt/scratch/local/blast/ncbi/nt":
            cmd_blast_db = " ".join(["makeblastdb",
                                     "-in",
                                     OTU_DATABASE,
                                     "-dbtype",
                                     "nucl"])
            logger.info("%s make blastdb command", cmd_blast_db)
            pipe = subprocess.run(cmd_blast_db, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True)
            logger.info("BLAST makde database stdout: %s", pipe.stdout)
            logger.info("BLAST makde database: %s", pipe.stderr)
    # run the blast
    xml_out = "%s_vs_%s.xml" % (fa_out_QC.split("qc")[0],
                                db_name.split(".fa")[0])
    cmd_blastrun = " ".join(["blastn",
                             "-db",
                             OTU_DATABASE,
                             "-query",
                             fa_out_QC,
                             "-num_threads",
                             args.threads,
                             "-evalue",
                             str(args.evalue),
                             "-outfmt 5",
                             "-out",
                             xml_out])
    logger.info("%s make cmd_blastrun command", cmd_blastrun)
    pipe = subprocess.run(cmd_blastrun, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
    logger.info("BLAST stdout: %s", pipe.stdout)
    logger.info("BLAST stderr: %s", pipe.stderr)
    ####################################################################
    # parse xml
    out_file_name = "%s_V_%s.txt" %(fa_out_QC.split("qc")[0],
                                    db_name.split(".fa")[0])
    cmd_parse = make_parse_cmd(out_file_name, args.mismatches)
    run_parse_cmd(cmd_parse, logger)
    mismatches = int(args.mismatches)
    while os.path.getsize(out_file_name) < 1:
        logger.warning("no hits at %d mismatches", mismatches)
        os.remove(out_file_name)
        mismatches = mismatches + 1
        logger.warning("going to increasing by 1. Mismatches = %d" % mismatches)
        cmd_parse = make_parse_cmd(out_file_name, mismatches)
        run_parse_cmd(cmd_parse, logger)
        if mismatches == 30:
            break
    remove_list = [fq_out, fa_out_all, fa_out, xml_out]
    if args.cleanup == "yes":
        for unwanted in remove_list:
            try:
                os.remove(unwanted)
                logger.warning("deleting: %s", unwanted)
            except:
                logger.info("could not find %s", unwanted)
    logger.info("Pipeline complete: %s", time.asctime())
