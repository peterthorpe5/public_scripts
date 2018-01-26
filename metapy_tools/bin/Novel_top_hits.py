#!/usr/bin/env python

#
# (c) The James Hutton Institute 2018
# Author: Peter Thorpe, Leighton Pritchard

import os
import subprocess
import sys
import shutil
import errno
import logging
import logging.handlers
import multiprocessing
import os
import time
import traceback
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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

VERSION = """Pycits/ metapy classify \n"""
if "--version" in sys.argv:
    print(VERSION)
    sys.exit(1)

cwd = os.getcwd()

def get_args():
    parser = argparse.ArgumentParser(description="Pipeline: cluster " +
                                     "data for metabarcoding ",
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("Nove")[0]
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--threads", dest='threads',
                          action="store", default="4",
                          type=str,
                          help="number of threads")

    optional.add_argument("--min_cluster_size",
                          dest='min_cluster_size',
                          action="store",
                          type=int,
                          default=200,
                          help="min number of seq in file to cosider looking into")

    optional.add_argument("-d", "--OTU_DB", dest='OTU_DB',
                          action="store",
                          default="nt",
                          type=str,
                          help="database of seq of to compare against")

    optional.add_argument("-e", "--evalue", dest='evalue',
                          action="store",
                          default=1e-50,
                          type=float,
                          help="evalue to filter results with")

    optional.add_argument("--num_seqs", dest='num_seqs',
                          action="store",
                          default="3",
                          type=str,
                          help="number of blast hits to get")

    optional.add_argument("-m", "--mismatches", dest='mismatches',
                          action="store",
                          default="20",
                          type=str,
                          help="number of mismatches to filter results with")

    optional.add_argument("--logfile",
                          dest="logfile",
                          action="store",
                          default="pipeline.log",
                          type=str,
                          help="Logfile name")
    optional.add_argument("--cleanup",
                          dest="cleanup",
                          action="store_true",
                          default="yes",
                          help="deletes most files the program creates: yes or no ")
    args = parser.parse_args()
    return args, file_directory


def count_seq_in_fasta(infasta):
    """func to count the seq in fasta file.
    Returns number of seq in file and first record in the
    file"""
    number_of_seq_count = 0
    for seq_record in SeqIO.parse(filename, "fasta"):
        if number_of_seq_count == 0:
            wanted_record = seq_record
        number_of_seq_count += 1
    return number_of_seq_count, wanted_record


def make_parse_cmd(outfile, mismatches):
    """function to return the parsing command"""
    cmd_parse = " ".join(["python",
                          os.path.join(FILE_DIRECTORY,
                                       "bin",
                                       "BLAST_xml_parser.py"),
                          "-i",
                          xml_out,
                          "-o",
                          outfile,
                          "-e",
                          str(args.evalue),
                          "-m",
                          str(mismatches)])
    return cmd_parse


def run_parse_cmd(cmd, logger):
    """func to run the parser command as this may be called a few times"""
    logger.info("%s make cmd_parse command", cmd)
    #  removed check-True
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)


def make_parse_cmd(outfile, mismatches):
    """function to return the parsing command"""
    cmd_parse = " ".join(["python",
                          os.path.join("$HOME",
                                       "public_scripts",
                                       "metapy_tools",
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


if __name__ == '__main__':
    # Set up logging
    args, FILE_DIRECTORY = get_args()
    PREFIX = "Novel_blast"
    logger = logging.getLogger('Novel_top_hit.py: %s' % time.asctime())
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
    OTU_DATABASE = args.OTU_DB
    logger.info("Command-line: %s", ' '.join(sys.argv))
    logger.info("Starting testing: %s", time.asctime())
    logger.info("using database: %s", OTU_DATABASE)
    if args.OTU_DB == "nt":
        logger.info("we are using nt. Tell the server where to find this")
        logger.info("you may need to alter this!")
        blast_export = "export BLASTDB=/mnt/scratch/local/blast/ncbi/"
        logger.info(blast_export)
        run_parse_cmd(blast_export, logger)
    else:
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
    # go through all fasta files in the directory
    for filename in os.listdir(".") :
        if not filename.endswith(".fasta"):
            continue
        # if the number of seq is greater than what we want
        number_of_seq_count, wanted_record = count_seq_in_fasta(filename)
        if number_of_seq_count > args.min_cluster_size:
            # write out a representative
            f_temp_out = open("temp.fasta", 'w')
            SeqIO.write(wanted_record, f_temp_out, "fasta")
            f_temp_out.close()
            # blast this representative against NR, or db of your choice.
            xml_out = filename.split(".fa")[0] + ".xml"
            cmd_blastrun = " ".join(["blastn",
                                     "-db",
                                     "/mnt/scratch/local/blast/ncbi/" + OTU_DATABASE,
                                     "-query",
                                     "temp.fasta",
                                     "-evalue",
                                     str(args.evalue),
                                     "-outfmt 5",
                                     "-max_target_seqs",
                                     args.num_seqs,
                                     "-num_threads",
                                     args.threads,
                                     "-out",
                                     xml_out])
            logger.info("%s make cmd_blastrun command", cmd_blastrun)
            pipe = subprocess.run(cmd_blastrun, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True)
            if pipe.stdout != "":
                logger.info("BLAST stdout: %s", pipe.stdout)
            if pipe.stderr != "":
                logger.info("BLAST stderr: %s", pipe.stderr)
            out_file_name = "%s_V_%s.txt" %(filename.split(".fa")[0],
                                            OTU_DATABASE)

            cmd_parse = make_parse_cmd(out_file_name, str(args.mismatches))
            run_parse_cmd(cmd_parse, logger)
            remove_list = [xml_out, "temp.fasta"]
            cleanup_option = args.cleanup.upper()
            if cleanup_option == "YES":
                for unwanted in remove_list:
                    try:
                        os.remove(unwanted)
                        logger.info("deleting: %s", unwanted)
                    except:
                        logger.info("could not find %s", unwanted)

            logger.info("finished %s" % filename)


