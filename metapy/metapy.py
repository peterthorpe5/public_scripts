#!/usr/bin/env python3
#
# metapy.py
#
# Code for script to identify OTUs from metabarcoding reads.
# runs multiple clustering programs
#
# (c) The James Hutton Institute 2016-2017
# Author: Leighton Pritchard, Peter Thorpe

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
import sklearn
import argparse
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from pycits.tools import convert_fq_to_fa, NotExecutableError, trim_seq,\
     dereplicate_name, check_OTU_db_abundance_val,\
     parse_tab_file_get_clusters, filter_sam_file, reformat_cdhit_clustrs,\
     reformat_sam_clusters, reformat_swarm_cls, reformat_blast6_clusters

from pycits.metapy_tools import decompress, compress,\
     last_exception, metapy_trim_seq, covert_chop_read, make_folder,\
     test_reads_exist_and_suffix, database_checker,\
     get_sizes, db_len_assembled_len_ok, stats_on_list_of_sizes,\
     plot_seq_len_histograms, full_illegal_charac_check

from pycits.Rand_index import pairwise_comparison_Rand

from pycits import tools, fastqc, trimmomatic, pear, error_correction,\
     flash, clean_up, swarm, seq_crumbs, bowtie_build, bowtie_map,\
     cd_hit, blast, vsearch, samtools_index, muscle

if sys.version_info[:1] != (3,):
    # e.g. sys.version_info(major=3, minor=6, micro=7,
    # releaselevel='final', serial=0)
    # break the program
    print ("currently using:", sys.version_info,
           "  version of python")
    raise ImportError("Python 3.x is now required for LTG.py")
    print ("did you activate the virtual environment?")
    print ("this is to deal with module imports")
    sys.exit(1)

VERSION = "Pycits classify OTU: v0.0.3"
if "--version" in sys.argv:
    print(VERSION)
    sys.exit(1)

#########################################################################


def get_args():
    parser = argparse.ArgumentParser(description="Pipeline: cluster " +
                                     "data for metabarcoding ",
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("metapy")[0]
    if not os.path.isfile(os.path.join(file_directory, "metapy.py")):
        file_directory = os.path.join(file_directory, "metapy")
    if not os.path.isfile(os.path.join(file_directory, "metapy.py")):
        print("Cannot locate correct path to metapy.py")
        sys.exit(1)                          
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--thread", dest='threads',
                          action="store", default="4",
                          type=str,
                          help="number of threads")

    # TODO: Why are these files being specified? They are not general.
    # for now they are here for testing. Will removed later
    # DESIGN: Attempting to collate everything here into a single cmd-line
    #         interface is a bad idea. We need to restructure this.
    optional.add_argument("-l", "--left", dest='left',
                          action="store",
                          default=os.path.join(file_directory,
                                               "data",
                                               "reads",
                                               "DNAMIX_S95_" +
                                               "L001_R1_001.fastq.gz"),
                          type=str,
                          help="left illumina reads, " +
                          "default is for tests")

    optional.add_argument("-r", "--right", dest='right',
                          action="store",
                          default=os.path.join(file_directory,
                                               "data",
                                               "reads",
                                               "DNAMIX_S95_" +
                                               "L001_R2_001.fastq.gz"),
                          type=str,
                          help="right illumina reads")

    optional.add_argument("-d", "--OTU_DB", dest='OTU_DB',
                          action="store",
                          default=os.path.join(file_directory,
                                               "data",
                                               "Phytophora_ITS_" +
                                               "database_v0.005.fasta"),
                          type=str,
                          help="right illumina reads")

    optional.add_argument("-a", "--assemble", dest='assemble',
                          action="store",
                          choices=["pear", "flash"],
                          help="program to assemble with " +
                          "flash or pear " +
                          "default: %(default)s", default="pear",
                          type=str)

    optional.add_argument("--adaptors", dest='adaptors',
                          action="store",
                          default=os.path.join(file_directory,
                                               "adapters",
                                               "TruSeq3-PE.fa"),
                          type=str,
                          help="adaptors for trimming. Can supply custom " +
                          " file if desired")

    optional.add_argument("--left_trim", dest='left_trim',
                          action="store",
                          default=53,
                          type=int,
                          help="left_trim for primers or conserved " +
                          "regions. Default 53 ")

    optional.add_argument("--right_trim", dest='right_trim',
                          action="store",
                          default=20,
                          type=int,
                          help="right_trim for primers or conserved " +
                          "regions.  Default 20 (5.8 region)")

    optional.add_argument("--phred", dest='phred',
                          action="store",
                          default="phred33",
                          type=str,
                          help="phred33 is default. " +
                          "Dont change unless sure")

    optional.add_argument("--cdhit_threshold", dest='cdhit_threshold',
                          action="store",
                          default="0.99",
                          type=str,
                          help="percentage identify for cd-hit " +
                          "Default -0.99")

    optional.add_argument("--swarm_d_value", dest='swarm_d_value',
                          action="store",
                          default=1,
                          type=int,
                          help="the difference d value for clustering " +
                          "in swarm. Default 1")

    optional.add_argument("--blastclust_threshold",
                          dest='blastclust_threshold',
                          action="store",
                          default=0.90,
                          type=float,
                          help="the threshold for blastclust clustering " +
                          " Default -S 0.90")

    optional.add_argument("--vesearch_threshold",
                          dest='vesearch_threshold',
                          action="store",
                          default=0.99,
                          type=float,
                          help="the threshold for vsearch clustering " +
                          " Default 0.99")

    optional.add_argument("--verbose", dest="verbose",
                          action="store_true",
                          default=False,
                          help="Report verbose output")

    optional.add_argument("-e", "--error_correction",
                          dest="Error_correction",
                          action="store_true",
                          default=False,
                          help="to perform Illumina error correction ")

    optional.add_argument("--align", dest="align",
                          action="store_true",
                          default=False,
                          help="to align clusters in the output " +
                          "you must have muscle in your PATH as muscle")

    optional.add_argument("--percent_identity",
                          dest="percent_identity",
                          action="store_true",
                          default=False,
                          help="blast the cluster to return " +
                          "pairwise percentage identity")

    optional.add_argument("--min_novel_cluster_threshold",
                          dest="min_novel_cluster_threshold",
                          type=str,
                          default="2",
                          help="min size of a cluster to consider as real " +
                          "anything smaller than this is ignored")

    optional.add_argument("--blastclust", dest="blastclust",
                          action="store", default="blastclust",
                          help="Path to blastclust ... If version already" +
                          "in PATH then leave blank")

    optional.add_argument("--muscle", dest="muscle",
                          action="store", default="muscle",
                          help="Path to MUSCLE... If version already" +
                          "in PATH then leave blank")

    optional.add_argument("--flash", dest="flash",
                          action="store", default="flash",
                          help="Path to flash... If version already" +
                          "in PATH then leave blank")

    optional.add_argument("--pear", dest="pear",
                          action="store", default="pear",
                          help="Path to pear... If version already" +
                          "in PATH then leave blank")

    optional.add_argument("--cd-hit-est", dest="cd_hit",
                          action="store", default="cd-hit-est",
                          help="Path to cd-hit-est... If version already" +
                          "in PATH then leave blank")

    optional.add_argument("--bowtie2", dest="bowtie2",
                          action="store", default="bowtie2",
                          help="Path to bowtie2... If version already" +
                          "in PATH then leave blank")

    optional.add_argument("--fastqc", dest="fastqc",
                          action="store", default="fastqc",
                          help="Path to fastqc... If version already" +
                          "in PATH then leave blank")

    optional.add_argument("--spades", dest="spades",
                          action="store", default="spades.py",
                          help="Path to spades.py... If version already" +
                          "in PATH then leave blank")

    optional.add_argument("--vsearch", dest="vsearch",
                          action="store", default="vsearch",
                          help="Path to vsearch... If version already" +
                          "in PATH then leave blank")

    optional.add_argument("--trimmomatic", dest="trimmomatic",
                          action="store", default="trimmomatic",
                          help="Path to trimmomatic... If version already" +
                          "in PATH then leave blank")

    optional.add_argument("--swarm", dest="swarm",
                          action="store", default="swarm",
                          help="Path to swarm... If version already" +
                          "in PATH then leave blank")

    optional.add_argument("--samtools", dest="samtools",
                          action="store", default="samtools",
                          help="Path to samtools... If version already" +
                          "in PATH then leave blank")

    optional.add_argument("--logfile",
                          dest="logfile",
                          action="store",
                          default="pipeline.log",
                          type=str,
                          help="Logfile name")

    optional.add_argument("--Run_blastclust",
                          dest="Run_blastclust",
                          action="store",
                          default="no",
                          type=str,
                          help="Run_blastclust yes or no, This is slow" +
                          " default no")

    optional.add_argument("--Run_cdhit",
                          dest="Run_cdhit",
                          action="store",
                          default="no",
                          type=str,
                          help="Run_cdhit yes or no, This is slow" +
                          " default no")

    optional.add_argument("--Run_dada2",
                          dest="Run_dada2",
                          action="store",
                          default="no",
                          type=str,
                          help="Run_dada2 yes or no, " +
                          " default no")

    optional.add_argument("--plot",
                          dest="plots",
                          action="store",
                          default="NO",
                          type=str,
                          help="create graphs of clusters yes or no, " +
                          "This is slow. Default no. ")

    optional.add_argument("--pvalue",
                          dest="pvalue",
                          action="store",
                          default=0.0000001,
                          type=float,
                          help="pvalue for comparing the database " +
                          "lengths versus the assembled lengths. At " +
                          "which point are the different?")

    optional.add_argument("--standard_deviation",
                          dest="std",
                          action="store",
                          default=2,
                          type=int,
                          help="standard_deviation threshold for " +
                          "comparing the assembled size versus the " +
                          "database sequence sizes. This is to check " +
                          "database is sensible for the data input")

    optional.add_argument("--cleanup",
                          dest="cleanup",
                          action="store",
                          type=str,
                          default="YES",
                          help="deletes most files the program creates " +
                          " yes/ no")
    optional.add_argument("--qc",
                          dest="qc",
                          action="store_true",
                          default=True,
                          help="performs QC at various stages. " +
                          "Turn off by: --qc False")

    optional.add_argument("-h", "--help",
                          action="help",
                          default=argparse.SUPPRESS,
                          help="Displays this help message"
                          " type --version for version")
    optional.add_argument('-v', '--version',
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
# test read set are defaults
LEFT_READS = test_reads_exist_and_suffix(args.left)
RIGHT_READS = test_reads_exist_and_suffix(args.right)
assert(LEFT_READS != RIGHT_READS)
READ_PREFIX = os.path.split(LEFT_READS)[-1].split("_R")[0]
PHREDSCORE = args.phred
PREFIX = READ_PREFIX
ADAPTERS = args.adaptors
OUTDIR_TRIM = "trimmed_reads"
OUTFILES = [os.path.join(OUTDIR_TRIM, PREFIX + suffix) for suffix in
            ("_paired_R1.fq.gz", "_unpaired_R1.fq.gz",
             "_paired_R2.fq.gz", "_unpaired_R2.fq.gz")]

OTU_DATABASE = args.OTU_DB
WORKING_DIR = os.getcwd()
# set this to false for now to not run it
SEQ_CRUMBS = False
CDHIT_THRESHOLD = str(args.cdhit_threshold)
if float(CDHIT_THRESHOLD) < 0.8:
    sys.exit("\nCDHIT_THRESHOLD must be less than 0.8\n")
SWARM_D_VALUE = args.swarm_d_value
if SWARM_D_VALUE > 2:
    WARNINGS = "swarm threshold is very 'loose' = %d\n" % SWARM_D_VALUE
VSEARCH_THRESHOLD = args.vesearch_threshold
ERROR_CORRECTION = args.Error_correction
ASSEMBLE_PROG = args.assemble

# left_trim 53 - this should be default??
LEFT_TRIM = args.left_trim
RIGHT_TRIM = args.right_trim
assert(READ_PREFIX == os.path.split(RIGHT_READS)[-1].split("_R")[0])
RESULT_FOLDER = PREFIX + "_RESULTS"
make_folder(RESULT_FOLDER, WORKING_DIR)
CLUSTER_FILES_FOR_RAND_INDEX = []
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
        flash.Flash(args.flash)
        tools_list.append("flash")
    except ValueError:
        Warning_out = Warning_out + "Flash not in path"
    try:
        error_correction.Error_Correction(args.spades)
        tools_list.append("error_correction")
    except ValueError:
        Warning_out = Warning_out + "spades.py not in path\n"
    try:
        vsearch.Vsearch(args.vsearch)
        tools_list.append("vsearch")
    except ValueError:
        Warning_out = Warning_out + "vsearch not in path\n"
    try:
        trimmomatic.Trimmomatic(args.trimmomatic)
        tools_list.append("trimmomatic")
    except ValueError:
        Warning_out = Warning_out + "trimmomatic not in path\n"
    try:
        swarm.Swarm(args.swarm)
        tools_list.append("swarm")
    except ValueError:
        Warning_out = Warning_out + "swarm not in path\n"
    try:
        samtools_index.Samtools_Index(args.samtools)
        tools_list.append("samtools")
    except ValueError:
        Warning_out = Warning_out + "samtools not in path\n"
    try:
        pear.Pear(args.pear)
        tools_list.append("pear")
    except ValueError:
        Warning_out = Warning_out + "pear not in path\n"
    try:
        muscle.Muscle(args.muscle)
        tools_list.append("muscle")
    except ValueError:
        Warning_out = Warning_out + "muscle not in path\n"
    try:
        fastqc.FastQC(args.fastqc)
        tools_list.append("fastqc")
    except ValueError:
        Warning_out = Warning_out + "fastqc not in path\n"
    try:
        cd_hit.Cd_hit(args.cd_hit)
        tools_list.append("cd-hit-est")
    except ValueError:
        Warning_out = Warning_out + "cd-hit-est not in path\n"
    try:
        bowtie_map.Bowtie2_Map(args.bowtie2)
        tools_list.append("bowtie2")
    except ValueError:
        Warning_out = Warning_out + "bowtie2 not in path\n"
    try:
        blast.Blastclust(args.blastclust)
        tools_list.append("blastclust")
    except ValueError:
        Warning_out = Warning_out + "blastclust not in path\n"
    return tools_list, Warning_out


# DESIGN: Why are you ignoring all the wrappers? This is hugely
#         counterproductive!
# DESIGN: If you write a Python script, and need to call code from another
#         script, that should be telling you that you need to write
#         a function in a module, and call it from both scripts. Using
#         a shell call to Python in your Python script is *very bad*!

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
    db_warnings = full_illegal_charac_check(OTU_DATABASE)
    db_error = False
    if db_warnings != "":
        logger.warning("DBwarnings: %s ", db_warnings)
        db_error = True
    # set up the DB and check if this is formatted correctly
    # function returns "ok", "ok" if passed.
    # else it returns:"Duplicate names found", seq_record
    # or "Duplicate sequence found", seq_record
    # or "Ill formatted fasta file", seq_record
    # depending on the problem
    value1, value2 = database_checker(OTU_DATABASE)
    if args.qc:
        if value1 != "ok":
            logger.warning("DATABASE CHECK FAILED.")
            logger.warning("Problem was: %s with %s" % (value1,
                                                        value2))
            logger.warning("check you database file")
            db_error = True
            # os._exit(0)  # should be kill here?
    # Get a list of tools in path!
    logger.info("checking which programs are in PATH")
    tools_list, Warning_out = check_tools_exist(WARNINGS)
    logger.info(Warning_out)

    ####################################################################
    # fastqc QC of raw reads
    if "fastqc" in tools_list:
        FASTQC_FOLDER = make_folder(os.path.join(RESULT_FOLDER,
                                                 PREFIX + "_fastqc"),
                                    WORKING_DIR)
        logger.info("starting fastqc")
        logger.info("made folder  %s",  FASTQC_FOLDER)
        qc = fastqc.FastQC(args.fastqc)
        qc_results = qc.run(LEFT_READS, PREFIX + FASTQC_FOLDER)
        logger.info("fastqc output: %s", qc_results.command)
        logger.info("fastqc stderr: %s", qc_results.stderr)

    ####################################################################
    # trimmomatic trm reads
    if "trimmomatic" in tools_list:
        TRIM_FOLDER = make_folder(os.path.join(RESULT_FOLDER,
                                               PREFIX + "_timmomatic"),
                                  WORKING_DIR)
        logger.info("starting trimmomatic testing")
        logger.info("made folder %s", TRIM_FOLDER)
        trim = trimmomatic.Trimmomatic(args.trimmomatic)

        logger.info("Trim reads by quality")
        parameters = trimmomatic.Parameters(threads=4)
        steps = trimmomatic.Steps(ILLUMINACLIP="{0}:2:30:10".format(ADAPTERS))
        results = trim.run(LEFT_READS, RIGHT_READS,
                           TRIM_FOLDER, PREFIX,
                           PHREDSCORE, parameters,
                           steps)
        logger.info("Trimming returned:", results)
        logger.info("Trimming command: %s", results.command)
        # get these exact file names from the named tuple
        logger.info("Trimming output: %s", results.stderr)
        LEFT_TRIMMED = results.outfileR1paired
        RIGHT_TRIMMED = results.outfileR2paired

    ####################################################################
    # error correction testing
    if "error_correction" in tools_list:
        if ERROR_CORRECTION:
            # we will error correct the read and reasign LEFT_TRIMMED
            # with the EC reads
            EC_FOLDER = make_folder(os.path.join(RESULT_FOLDER,
                                                 PREFIX + "_EC"),
                                    WORKING_DIR)
            error_corr = error_correction.Error_Correction(args.spades)
            logger.info("error correction using Bayes hammer")
            logger.info("made folder %s", EC_FOLDER)
            EC_results = error_corr.run(LEFT_TRIMMED,
                                        RIGHT_TRIMMED,
                                        THREADS,
                                        EC_FOLDER)
            # get these exact file names from the named tuple
            logger.info("error correc output: %s", EC_results.stderr)
            # assign the EC to the LEFT and RIGHT trimmed variable
            LEFT_TRIMMED = EC_results.Left_read_correct
            RIGHT_TRIMMED = EC_results.right_read_correct
            logger.info("Trimmed reads have been error corrected " +
                        "using Bayes hammer")

    ####################################################################
    # PEAR testing - assemble
    if "pear" in tools_list and ASSEMBLE_PROG == "pear":
        logger.info("assembly using PEAR")
        if ERROR_CORRECTION:
            suffix = "_PEAR_EC"
            logger.info("PEAR will use trimmed and error corrected reads")
        else:
            suffix = "_PEAR"
        ASSEMBLY_FOLDER = make_folder(os.path.join(RESULT_FOLDER,
                                                   PREFIX + suffix),
                                      WORKING_DIR)
        # call the class
        assemble = pear.Pear(args.pear)
        results_pear = assemble.run(LEFT_TRIMMED,
                                    RIGHT_TRIMMED,
                                    THREADS,
                                    ASSEMBLY_FOLDER,
                                    PREFIX)
        logger.info("PEAR returned:", results_pear)
        logger.info("PEAR command: %s", results_pear.command)
        logger.info("PEAR output: %s", results_pear.stderr)
        # this is from the named tuple
        ASSEMBLED = results_pear.outfileassembled
        # call the function
        logger.info("chopping the reads %s", ASSEMBLED)
        covert_chop_read(ASSEMBLED, LEFT_TRIM, RIGHT_TRIM)

    ####################################################################
    # FLASH testing - assemble
    if "flash" in tools_list and ASSEMBLE_PROG == "flash":
        logger.info("assembly using FLASH")
        if ERROR_CORRECTION:
            suffix = "_FLASH_EC"
            logger.info("FLASH will use trimmed and error corrected reads")
        else:
            suffix = "_FLASH"
        ASSEMBLY_FOLDER = make_folder(os.path.join(RESULT_FOLDER,
                                                   PREFIX + suffix),
                                      WORKING_DIR)
        assemble = flash.Flash(args.flash)
        logger.info("assembly using Flash")
        results_flash = assemble.run(LEFT_TRIMMED,
                                     RIGHT_TRIMMED,
                                     THREADS,
                                     ASSEMBLY_FOLDER,
                                     PREFIX)
        logger.info("Flash returned:", results_flash)
        logger.info("Flash command: %s", results_flash.command)
        logger.info("Flash output: %s", results_flash.stderr)
        ASSEMBLED = results_flash.outfileextended
        # call the function from tools
        logger.info("chopping the reads %s", ASSEMBLED)
        covert_chop_read(ASSEMBLED, LEFT_TRIM, RIGHT_TRIM)

    ####################################################################
    # convert format using seq_crumbs
    if SEQ_CRUMBS:
        format_change = seq_crumbs.Convert_Format("convert_format",
                                                  logger)
        format_change.run(ASSEMBLED,
                          ASSEMBLY_FOLDER,
                          logger)
    ####################################################################
    # check the assembled sizes are sensible against the db sizes.
    # are the db sequences significantly longer than the assembled sizes?
    # (db_fa, assembled_fa, sd=3)
    if args.qc:
        logger.info("QC: Checking your assembled seq sizes against db")
        # call fucntion from metapy_tools to get seq lens as list
        db_lens = get_sizes(OTU_DATABASE)
        assemb_lens = get_sizes(ASSEMBLED + ".bio.chopp" + "ed.fasta")
        # plot the seq len distributions
        plot_seq_len_histograms(RESULT_FOLDER, db_lens, assemb_lens)
        # call stats function to compare these:
        stats_data = stats_on_list_of_sizes(db_lens, assemb_lens)
        # metabarcoding will always have a skew. Especially if one
        # species is highly or only present.
        as_skew, db_skew, ttest, Man_u_value,\
                 Man_p_value = stats_data.split("\t")
        skew_t = "\t".join(["assembled_skew: %s" % as_skew,
                            "database_skew: %s" % db_skew,
                            "Mann_whitney U test: %s" % Man_p_value])
        logger.info("Assembled seq lengths versus database lengths")
        logger.info("%s", skew_t)
        # call the function to perform a simple mean versus
        # number of std test to find the problem
        size_test, error = db_len_assembled_len_ok(db_lens,
                                                   assemb_lens,
                                                   args.std)
        errtype, db_mean, db_sd, assemb_mean,\
                     assemb_sd = error.split("\t")
        dbstats = skew_t + "\tdb_mean= %s\tdb_stdev= %s\t" % (str(db_mean),
                                                              str(db_sd))
        stats = "assem_mean = %s , assem_stdev = %s" % (str(assemb_mean),
                                                        str(assemb_sd))
        if float(Man_p_value) < args.pvalue:
            if size_test == "fail":
                terminate = " ".join(["The assembled size of your reads is",
                                      "significantly different to your",
                                      "database. You need to adjust your",
                                      "DB sequences to that of the region",
                                      "you sequenced. \n"])
                if errtype == "-":
                    shorter = " ".join(["your assembled sequences are",
                                        "significantly shorter than db\n"])
                    error_out = "%s%s%s%s" % (shorter, terminate,
                                              dbstats, stats)
                    logger.warning("%s", error_out)
                    # KILL the program?
                    # sys.exit(error_out)
                else:
                    longer = " ".join(["your assemebled sequences are",
                                      "significantly longer than db\n"])
                    error_out = "%s%s%s%s" % (longer, terminate,
                                              dbstats, stats)
                    logger.warning("Warning: %s", error_out)
                    # KILL the program?
                    # sys.exit(error_out)
        logger.info("QC passed on sequences: %s\t%s", dbstats, stats)
    #######################################################################
    # first cat the db and EC, trimmed reads.
    cat_cmd = ["cat", OTU_DATABASE,
               ASSEMBLED + ".bio.chopped.fasta",
               ">",
               "assembled_fa_and_OTU_db.fasta"]
    cat_cmd = ' '.join(cat_cmd)
    logger.info("combine these files %s", cat_cmd)
    pipe = subprocess.run(cat_cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)

    ####################################################################
    # deduplicate reads:
    # use the function in tools dereplicate_name()
    if "swarm" in tools_list:
        logger.info("deduplicating reads")
        # dereplicate_name (infasta, db_out, out_fasta)
        dereplicate_name(ASSEMBLED + ".bio.chopped.fasta",
                         "db_old_to_new_names.txt",
                         ASSEMBLED + "for_swarm.fasta")

    ####################################################################
    # SWARM testing - assemble
    if "swarm" in tools_list:
        if db_error:
            logger.warning("error was found in your database. " +
                           "Swarm will most likely fail here. \n" +
                           "Look at the log file to find error")
        swarm_folder_name = PREFIX + "_Swarm_d%d" % SWARM_D_VALUE
        swarm_parameters = swarm.Parameters(t=1, d=SWARM_D_VALUE)
        SWARM_FOLDER = make_folder(os.path.join(RESULT_FOLDER,
                                                swarm_folder_name),
                                   WORKING_DIR)
        cluster = swarm.Swarm(args.swarm)
        assembled_fa_reads = ASSEMBLED + "for_swarm.fasta"
        logger.info("clustering with Swarm")
        # need to check the OTU database has abundance value
        # call the function from tools
        logger.info("OTU was %s", OTU_DATABASE)
        OTU_DATABASE_SWARM = check_OTU_db_abundance_val(OTU_DATABASE)
        logger.info("OTU is %s", OTU_DATABASE_SWARM)
        # need to cat the assembled_fasta with the database
        cat_cmd = ["cat",
                   OTU_DATABASE_SWARM,
                   assembled_fa_reads,
                   ">",
                   "assembled_reads_and_OTU_db.fasta"]
        cat_cmd = ' '.join(cat_cmd)
        logger.info("combine these files %s" % cat_cmd)
        pipe = subprocess.run(cat_cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        cluster_outdata = cluster.run("assembled_reads_and_OTU_db.fasta",
                                      SWARM_FOLDER,
                                      swarm_parameters)
        # use named tuple to get the outpfile name
        SWARM_OUT = cluster_outdata.outfilename

        logger.info("swarm returned using error corrected reads: ",
                    cluster_outdata)
        logger.info("swarm command: %s", cluster_outdata.command)
        logger.info("swarm output: %s", cluster_outdata.stderr)

        ################################################################
        # recode cluster output
        logger.info("renaming swarm output")
        # parse_tab_file_get_clusters(filename1, database, out_file)
        parse_tab_file_get_clusters(SWARM_OUT,
                                    OTU_DATABASE_SWARM,
                                    "db_old_to_new_names.txt",
                                    SWARM_OUT + "RENAMED_abundance",
                                    True)
        parse_tab_file_get_clusters(SWARM_OUT,
                                    OTU_DATABASE,
                                    "db_old_to_new_names.txt",
                                    SWARM_OUT + "RENAMED")
        # reformat_swarm_cls(swarm, db, db_and_reads, outfile)
        logger.info("reformatting swarm output for post analysis")

        reformat_swarm_cls(SWARM_OUT + "RENAMED",
                           OTU_DATABASE,
                           "assembled_fa_and_OTU_db.fasta",
                           SWARM_OUT + "for_R",
                           False)
        # add this file for Rand index comparison later
        CLUSTER_FILES_FOR_RAND_INDEX.append(SWARM_OUT + "for_R")
        cmd_s = ["python",
                 os.path.join(FILE_DIRECTORY,
                              "post_analysis",
                              "get_results_from_cluster_and_novel_" +
                              "clusterings.py"),
                 "-f", ASSEMBLED + ".bio.chopped.fasta",
                 "--all_fasta", "assembled_reads_and_OTU_db.fasta",
                 "--seq_db", OTU_DATABASE_SWARM,
                 "--min_novel_cluster_threshold",
                 args.min_novel_cluster_threshold,
                 "--left", LEFT_READS,
                 "--right", RIGHT_READS,
                 "--Name_of_project",
                 os.path.join(SWARM_FOLDER, "clusters"),
                 "--in",
                 SWARM_OUT + "RENAMED_abundance",
                 "--difference", str(SWARM_D_VALUE),
                 "-o",
                 os.path.join(SWARM_FOLDER,
                              "%s_swarm_results_%d.RESULTS" %
                              (PREFIX, SWARM_D_VALUE)),
                 "--old_to_new", "db_old_to_new_names.txt"]
        swarm_result = os.path.join(SWARM_FOLDER,
                                    "%s_swarm_results_%d.RESULTS" %
                                    (PREFIX, SWARM_D_VALUE))
        RESULTS.append("swarm\t%s" % swarm_result)
        cmd_s = ' '.join(cmd_s)
        if args.align:
            logger.info("going to align the cluster. Will take ages!")
            cmd_s = cmd_s + " --align True"
        if args.percent_identity:
            cmd_s = cmd_s + " --blast True"
        logger.info("%s = post analysis command", cmd_s)
        pipe = subprocess.run(cmd_s, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        if args.plots.upper() == "YES":
            logger.info("graphically represent swarm clusters")
            plot_cmd = ["python",
                        os.path.join(FILE_DIRECTORY,
                                     "bin",
                                     "draw_bar_chart_of_clusters.py"),
                        "-i",
                        SWARM_OUT + "RENAMED_abundance"
                        " --db",
                        OTU_DATABASE_SWARM]
            plot_cmd = ' '.join(plot_cmd)
            logger.info("plotting command = %s", plot_cmd)
            pipe = subprocess.run(plot_cmd, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True)

    #####################################################################
    # run cd hit
    # first cat the db and EC (if done with ec), trimmed reads.
    cat_cmd = ["cat", OTU_DATABASE,
               ASSEMBLED + ".bio.chopped.fasta",
               ">",
               "assembled_fa_and_OTU_db.fasta"]
    cat_cmd = ' '.join(cat_cmd)
    logger.info("combine these files: %s", cat_cmd)
    pipe = subprocess.run(cat_cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
##    if "cd-hit-est" in tools_list:
##        # default dont run cd hit!
##        if args.Run_cdhit.upper == "YES":
##            cd_folder_name = PREFIX + "_cd_hit_%s" % CDHIT_THRESHOLD
##            CDHIT_FOLDER = make_folder(os.path.join(RESULT_FOLDER,
##                                                    cd_folder_name),
##                                       WORKING_DIR)
##            cluster = cd_hit.Cd_hit(args.cd_hit)
##            results = cluster.run("assembled_fa_and_OTU_db.fasta",
##                                  THREADS,
##                                  CDHIT_THRESHOLD,
##                                  CDHIT_FOLDER,
##                                  PREFIX)
##            logger.info("cdhit: %s", results.command)
##            logger.info("reformatting cd hit output")
##            reformat_cdhit_clustrs(results.clusters,
##                                   results.clusters + "_1_line_per",
##                                   results.clusters + "for_R")
##            # add this file for Rand index comparison later
##            CLUSTER_FILES_FOR_RAND_INDEX.append(results.clusters + "for_R")
##            # analyse the clusters
##            cmd_c = ["python",
##                     os.path.join(FILE_DIRECTORY,
##                                  "post_analysis",
##                                  "get_results_from_cluster_and_novel_" +
##                                  "clusterings_cd_hit.py"),
##                     "-f", ASSEMBLED + ".bio.chopped.fasta",
##                     "--all_fasta", "assembled_fa_and_OTU_db.fasta",
##                     "--seq_db", OTU_DATABASE,
##                     "--min_novel_cluster_threshold",
##                     args.min_novel_cluster_threshold,
##                     "--left", LEFT_READS,
##                     "--right", RIGHT_READS,
##                     "--Name_of_project",
##                     os.path.join(CDHIT_FOLDER, "clusters"),
##                     "--in",
##                     results.clusters + "_1_line_per",
##                     "--difference", CDHIT_THRESHOLD,
##                     "-o",
##                     os.path.join(CDHIT_FOLDER,
##                                  "%s_cdhit_results_%s.RESULTS" %
##                                  (PREFIX, str(CDHIT_THRESHOLD))),
##                     "--old_to_new", "db_old_to_new_names.txt"]
##
##            cd_hit_result = os.path.join(CDHIT_FOLDER,
##                                         "%s_cdhit_results_%s.RESULTS" %
##                                         (PREFIX, str(CDHIT_THRESHOLD)))
##            RESULTS.append("cdhit\t%s" % cd_hit_result)
##
##            cmd_c = ' '.join(cmd_c)
##            if args.align:
##                logger.info("going to align the cluster. Will take ages!")
##                cmd_c = cmd_c + " --align True"
##            if args.percent_identity:
##                cmd_c = cmd_c + " --blast True"
##            logger.info("%s = post analysis command", cmd_c)
##            pipe = subprocess.run(cmd_c, shell=True,
##                                  stdout=subprocess.PIPE,
##                                  stderr=subprocess.PIPE,
##                                  check=True)
##            if args.plots.upper() == "YES":
##                logger.info("graphically represent cdhit clusters")
##                plot_cmd = ["python",
##                            os.path.join(FILE_DIRECTORY,
##                                         "bin",
##                                         "draw_bar_chart_of_clusters.py"),
##                            "-i",
##                            results.clusters + "_1_line_per",
##                            " --db",
##                            OTU_DATABASE]
##                plot_cmd = ' '.join(plot_cmd)
##                logger.info("plotting command = %s", plot_cmd)
##                pipe = subprocess.run(plot_cmd, shell=True,
##                                      stdout=subprocess.PIPE,
##                                      stderr=subprocess.PIPE,
##                                  check=True)

    ####################################################################
    # run vsearch
    if "vsearch" in tools_list:
        v_n = PREFIX + "_Vsearch_us_global_%.2f" % VSEARCH_THRESHOLD
        VSEARCH_FOLDER = make_folder(os.path.join(RESULT_FOLDER,
                                                  v_n),
                                     WORKING_DIR)
        OUTFILE_DEREP = os.path.join(VSEARCH_FOLDER,
                                     "vsearch_derep.fasta")
        OUTFILE_CLUSTER_UC = os.path.join(VSEARCH_FOLDER,
                                          PREFIX + "vsearch_cluster.uc")
        OUTFILE_CLUSTER_B6 = os.path.join(VSEARCH_FOLDER,
                                          PREFIX +
                                          "vsearch_cluster.blast6")

        # PARAMETERS
        CLUSTER_PARAMS = {'--blast6out': OUTFILE_CLUSTER_B6,
                          '--id': VSEARCH_THRESHOLD,
                          '--db': OTU_DATABASE,
                          '--threads': THREADS}

        vsearch_exe = vsearch.Vsearch(args.vsearch)
        Result_derep = vsearch_exe.run('--derep_fulllength',
                                       ASSEMBLED + ".bio.chopped.fasta",
                                       ASSEMBLED + "drep.vsearch.fasta")
        logger.info("vsearch derep: %s", Result_derep.command)
        logger.info("vsearch clustering")
        v_cluster = vsearch_exe.run('--usearch_global',
                                    ASSEMBLED + ".bio.chopped.fasta",
                                    OUTFILE_CLUSTER_UC,
                                    CLUSTER_PARAMS)
        logger.info("vsearch cluster: %s", v_cluster.command)
        reformat_blast6_clusters(v_cluster.outfile_b6,
                                 "assembled_fa_and_OTU_db.fasta",
                                 v_cluster.outfile_b6 + "for_R")
        # add this file for Rand index comparison later
        CLUSTER_FILES_FOR_RAND_INDEX.append(v_cluster.outfile_b6 + "for_R")

        cmd_v = ["python",
                 os.path.join(FILE_DIRECTORY,
                              "post_analysis",
                              "get_results_from_cluster_and_novel_" +
                              "clusterings_cd_hit.py"),
                 "-f", ASSEMBLED + ".bio.chopped.fasta",
                 "--all_fasta", "assembled_fa_and_OTU_db.fasta",
                 "--seq_db", OTU_DATABASE,
                 "--min_novel_cluster_threshold",
                 args.min_novel_cluster_threshold,
                 "--left", LEFT_READS,
                 "--right", RIGHT_READS,
                 "--Name_of_project",
                 os.path.join(VSEARCH_FOLDER, "clusters"),
                 "--in",
                 v_cluster.outfile_b6 + "for_R_1_line",
                 "--difference", str(VSEARCH_THRESHOLD),
                 "-o",
                 os.path.join(VSEARCH_FOLDER,
                              "%s_vsearch_results_%s.RESULTS" %
                              (PREFIX, str(VSEARCH_THRESHOLD))),
                 "--old_to_new", "db_old_to_new_names.txt"]

        v_result = os.path.join(VSEARCH_FOLDER,
                                "%s_vsearch_results_%s.RESULTS" %
                                (PREFIX, str(VSEARCH_THRESHOLD)))
        RESULTS.append("vsearch\t%s" % v_result)

        cmd_v = ' '.join(cmd_v)
        if args.align:
            logger.info("going to align the cluster. Will take ages!")
            cmd_v = cmd_v + " --align True"
        if args.percent_identity:
            cmd_v = cmd_v + " --blast True"
        logger.info("%s = post analysis command", cmd_v)
        pipe = subprocess.run(cmd_v, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        if args.plots.upper() == "YES":
            logger.info("graphically represent vsearch default clusters")
            plot_cmd = ["python",
                        os.path.join(FILE_DIRECTORY,
                                     "bin",
                                     "draw_bar_chart_of_clusters.py"),
                        "-i",
                        v_cluster.outfile_b6 + "for_R_1_line",
                        " --db",
                        OTU_DATABASE]
            plot_cmd = ' '.join(plot_cmd)
            logger.info("plotting command = %s", plot_cmd)
            pipe = subprocess.run(plot_cmd, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True)

    ###################################################################
    # deduplicate reads:
    # use the function in tools dereplicate_name()
    # run vsearch
    if "vsearch" in tools_list:
        # PARAMETERS
        vn2 = PREFIX + "_Vsearch_clustfast_%.2f" % VSEARCH_THRESHOLD
        VSEARCH_FOLDER = make_folder(os.path.join(RESULT_FOLDER,
                                                  vn2),
                                     WORKING_DIR)
        OUTFILE_CLUSTER_FAST_UC = os.path.join(VSEARCH_FOLDER,
                                               PREFIX +
                                               "vsearch_cluster_fast.uc")
        OUTFILE_CLUSTER_FAST_B6 = os.path.join(VSEARCH_FOLDER,
                                               PREFIX +
                                               "vsearch_cluster_fast." +
                                               "blast6")
        OUTFILE_CLUSTER_FAST_MSA = os.path.join(VSEARCH_FOLDER,
                                                PREFIX +
                                                "vsearch_cluster_fast" +
                                                "_msa.fasta")
        OUTFILE_CLUSTER_FAST_CENTROIDS = os.path.join(VSEARCH_FOLDER,
                                                      PREFIX +
                                                      "vsearch_cluster" +
                                                      "_fast_centroids." +
                                                      "fasta")
        OUTFILE_CLUSTER_FAST_CONSENSUS = os.path.join(VSEARCH_FOLDER,
                                                      PREFIX +
                                                      "vsearch_cluster" +
                                                      "_fast_consensus." +
                                                      "fasta")
        CLUSTER_FAST_PARAMS = {"--id": VSEARCH_THRESHOLD,
                               "--centroids": OUTFILE_CLUSTER_FAST_CENTROIDS,
                               "--msaout": OUTFILE_CLUSTER_FAST_MSA,
                               "--consout": OUTFILE_CLUSTER_FAST_CONSENSUS,
                               "--db": OTU_DATABASE,
                               "--threads": THREADS,
                               "--blast6out": OUTFILE_CLUSTER_FAST_B6}
        logger.info("deduplicating reads")
        # derep the database for vsearch
        db_derep = os.path.join(VSEARCH_FOLDER, PREFIX + "for_vsearch.fasta")
        derep_db = vsearch_exe.run('--derep_fulllength',
                                   OTU_DATABASE,
                                   db_derep)
        # derep the name, True means it is for vsearch
        dereplicate_name(ASSEMBLED + ".bio.chopped.fasta",
                         "db_old_to_new_names_vsearch.txt",
                         ASSEMBLED + "for_vsearch.fasta",
                         True)
        # cat the derep reads and derep db together
        cat_cmd = ["cat", derep_db.outfilename,
                   ASSEMBLED + "for_vsearch.fasta",
                   ">",
                   "assembled_fa_and_OTU_db_vesearch.fasta"]
        cat_cmd = ' '.join(cat_cmd)
        logger.info("merge these files %s", cat_cmd)
        pipe = subprocess.run(cat_cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        # use vsearch to get concensus and aligned clusters
        fast_clust = vsearch_exe.run("--cluster_fast",
                                     "assembled_fa_and_OTU_db_vesearch.fasta",
                                     OUTFILE_CLUSTER_FAST_UC,
                                     CLUSTER_FAST_PARAMS)

        reformat_blast6_clusters(fast_clust.outfile_b6,
                                 "assembled_fa_and_OTU_db.fasta",
                                 fast_clust.outfile_b6 + "for_R")
        logger.info("vsearch cluster: %s", fast_clust.command)
        CLUSTER_FILES_FOR_RAND_INDEX.append(fast_clust.outfile_b6 + "for_R")

        cmd_F = ["python",
                 os.path.join(FILE_DIRECTORY,
                              "post_analysis",
                              "get_results_from_cluster_and_novel_" +
                              "clusterings_cd_hit.py"),
                 "-f", ASSEMBLED + ".bio.chopped.fasta",
                 "--all_fasta", "assembled_fa_and_OTU_db.fasta",
                 "--seq_db", OTU_DATABASE,
                 "--min_novel_cluster_threshold",
                 args.min_novel_cluster_threshold,
                 "--left", LEFT_READS,
                 "--right", RIGHT_READS,
                 "--Name_of_project",
                 os.path.join(VSEARCH_FOLDER, "clusters_clusterfast"),
                 "--in",
                 fast_clust.outfile_b6 + "for_R_1_line",
                 "--difference", str(VSEARCH_THRESHOLD),
                 "-o",
                 os.path.join(VSEARCH_FOLDER,
                              "%s_vserach_clu_fasta_results_%s.RESULTS" %
                              (PREFIX, str(VSEARCH_THRESHOLD))),
                 "--old_to_new",
                 "db_old_to_new_names.txt"]

        v_f_result = os.path.join(VSEARCH_FOLDER,
                                  "%s_vserach_clu_fasta_results_%s.RESULTS" %
                                  (PREFIX, str(VSEARCH_THRESHOLD)))

        RESULTS.append("vse_faclus\t%s" % v_f_result)

        cmd_F = ' '.join(cmd_F)
        if args.align:
            logger.info("going to align the cluster. Will take ages!")
            cmd_F = cmd_F + " --align True"
        if args.percent_identity:
            cmd_F = cmd_F + " --blast True"
        logger.info("%s = post analysis command", cmd_v)
        pipe = subprocess.run(cmd_F, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        if args.plots.upper() == "YES":
            logger.info("graphically represent vsearch clu_fasta clusters")
            plot_cmd = ["python",
                        os.path.join(FILE_DIRECTORY,
                                     "bin",
                                     "draw_bar_chart_of_clusters.py"),
                        "-i",
                        fast_clust.outfile_b6 + "for_R_1_line",
                        " --db",
                        OTU_DATABASE]
            plot_cmd = ' '.join(plot_cmd)
            logger.info("plotting command = %s", plot_cmd)
            pipe = subprocess.run(plot_cmd, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True)

    #####################################################################
    # MAP THE READS WITH BOWTIE
    if "bowtie2" in tools_list:
        BOWTIE_FOLDER = make_folder(os.path.join(RESULT_FOLDER,
                                                 PREFIX + "_bowtie"),
                                    WORKING_DIR)
        index = bowtie_build.Bowtie2_Build("bowtie2-build")
        results = index.run(OTU_DATABASE, "OTU")
        logger.info("bowtie build %s", results.command)
        mapper = bowtie_map.Bowtie2_Map("bowtie2")
        otu = os.path.split(OTU_DATABASE)[-1]
        bowtie_out = os.path.join(BOWTIE_FOLDER,
                                  "_vs_".join([PREFIX,
                                               otu.split(".f")[0]]) + ".sam")
        # e.g.  result = bt2_map.run(READS, FA_INDEX,
        #  outfilename, THREADS, fasta=True)
        results = mapper.run(ASSEMBLED + ".bio.chopped.fasta",
                             "OTU",
                             bowtie_out,
                             THREADS,
                             fasta=True)
        logger.info("bowtie map %s", results.command)
        logger.info("pysam to filter the mapping")
        samfile = pysam.AlignmentFile(results.sam, "r")
        # call the function from tools
        # this return matches with zero mismatches, but not to be interpreted
        # as a perfect match?!?!
        cig_list, matches = filter_sam_file(results.sam,
                                            (os.path.join(BOWTIE_FOLDER,
                                                          "pysam_perfect_" +
                                                          "cigar.txt")))
        logger.info("pysam found %d perfect(?) cigar MATCHES" % (len(cig_list) -1))
        # using grep to get perfect matches:
        grep_cmd = ' '.join(['cat',
                             results.sam,
                             '|',
                             'grep',
                             '"AS:i:0"',
                             '>',
                             results.sam + "perfect_map"])
        logger.info("grep for perfect reads %s" % grep_cmd)
        pipe = subprocess.run(grep_cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        # sort these and make them unique.
        # for info only, this next command is not needed
        grep_cmd = ' '.join(['cat',
                             results.sam,
                             '|',
                             'grep',
                             '"AS:i:0"',
                             '|',
                             'cut -f3',
                             '|',
                             'sort'
                             '|',
                             'uniq'
                             '>',
                             results.sam + "perfect_map_name"])
        logger.info("grep for perfect reads %s" % grep_cmd)
        pipe = subprocess.run(grep_cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        perfect_count = 0
        with open(results.sam + "perfect_map_name") as perfect_names:
            data = perfect_names.read().split("\n")
            for name in data:
                perfect_count = perfect_count + 1
                logger.info("%s = perfect match", name)
        logger.info("%d = number perfect match", perfect_count)
        reformat_sam_clusters(results.sam + "perfect_map",
                              "assembled_fa_and_OTU_db.fasta",
                              results.sam +
                              "for_R")
        # add this file for Rand index comparison later
        CLUSTER_FILES_FOR_RAND_INDEX.append(results.sam + "for_R")

        cmd_b = ["python",
                 os.path.join(FILE_DIRECTORY,
                              "post_analysis",
                              "get_results_from_cluster_and_novel_" +
                              "clusterings_cd_hit.py"),
                 "-f", ASSEMBLED + ".bio.chopped.fasta",
                 "--all_fasta", "assembled_fa_and_OTU_db.fasta",
                 "--seq_db", OTU_DATABASE,
                 "--min_novel_cluster_threshold",
                 args.min_novel_cluster_threshold,
                 "--left", LEFT_READS,
                 "--right", RIGHT_READS,
                 "--Name_of_project",
                 os.path.join(BOWTIE_FOLDER, "clusters_perfect_map"),
                 "--in",
                 results.sam + "for_R_1_line",
                 "--difference", "perfect_map",
                 "-o",
                 os.path.join(BOWTIE_FOLDER,
                              PREFIX + "_bowtie_perfect.RESULTS"),
                 "--old_to_new", "db_old_to_new_names.txt"]
        bow_result = os.path.join(BOWTIE_FOLDER,
                                  PREFIX + "_bowtie_perfect.RESULTS")

        RESULTS.append("bowtie\t%s" % bow_result)

        cmd_b = ' '.join(cmd_b)
        logger.info("%s = post analysis command" % cmd_b)
        pipe = subprocess.run(cmd_b, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        if args.plots.upper() == "YES":
            logger.info("graphically represent bowtie clusters")
            plot_cmd = ["python",
                        os.path.join(FILE_DIRECTORY,
                                     "bin",
                                     "draw_bar_chart_of_clusters.py"),
                        "-i",
                        results.sam + "for_R_1_line",
                        " --db",
                        OTU_DATABASE]
            plot_cmd = ' '.join(plot_cmd)
            logger.info("plotting command = %s", plot_cmd)
            pipe = subprocess.run(plot_cmd, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True)

    #####################################################################
    #####################################################################
    # DADA2 pipeline - put sequences through Swarm
    if "swarm" in tools_list:
        # default not to run dada2
        if args.Run_dada2.upper() == "YES":
            logger.info("Running DADA2 pipeline")
            dada2_folder_name = PREFIX + "_DADA2_Swarm_d%d" % SWARM_D_VALUE
            swarm_parameters = swarm.Parameters(t=1, d=SWARM_D_VALUE)
            DADA2_FOLDER = make_folder(os.path.join(RESULT_FOLDER,
                                                    dada2_folder_name),
                                       WORKING_DIR)
            # run DADA2 pipeline by calling program from bin.
            # This will produce a fasta will final QC seq in it
            cmd_dada = ["python",
                        os.path.join(FILE_DIRECTORY,
                                     "bin",
                                     "wrap_DADA2.py"),
                        "-d", OTU_DATABASE_SWARM,
                        "--left", LEFT_READS,
                        "--right", RIGHT_READS]
            cmd_dada = ' '.join(cmd_dada)
            logger.info(cmd_dada)
            pipe = subprocess.run(cmd_dada, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True)

            # set up an instance of swarm for clustering later
            cluster = swarm.Swarm(args.swarm)
            # call the function from tools
            logger.info("OTU was %s", OTU_DATABASE)
            OTU_DATABASE_SWARM = check_OTU_db_abundance_val(OTU_DATABASE)
            logger.info("OTU is %s", OTU_DATABASE_SWARM)
            # need to cat the assembled_fasta with the database
            cat_cmd = ["cat",
                       OTU_DATABASE_SWARM,
                       READ_PREFIX + "_DADA2.fasta",
                       ">",
                       "dada2_seq_and_OTU_db.fasta"]
            cat_cmd = ' '.join(cat_cmd)
            logger.info("combine these files %s" % cat_cmd)
            pipe = subprocess.run(cat_cmd, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True)
            cluster_outdata = cluster.run("dada2_seq_and_OTU_db.fasta",
                                          DADA2_FOLDER,
                                          swarm_parameters)
            # use named tuple to get the outpfile name
            SWARM_OUT = cluster_outdata.outfilename

            logger.info("swarm command: %s", cluster_outdata.command)
            logger.info("swarm output: %s", cluster_outdata.stderr)

            ################################################################
            # recode cluster output
            logger.info("reformatting swarm output for post analysis")

            reformat_swarm_cls(SWARM_OUT,
                               OTU_DATABASE_SWARM,
                               "dada2_seq_and_OTU_db.fasta",
                               SWARM_OUT + "for_R",
                               False)
            # add this file for Rand index comparison later
            CLUSTER_FILES_FOR_RAND_INDEX.append(SWARM_OUT + "for_R")
            cmd_s = ["python",
                     os.path.join(FILE_DIRECTORY,
                                  "post_analysis",
                                  "get_results_from_cluster_and_novel_" +
                                  "clusterings_dada2.py"),
                     "-f", "dada2_seq_and_OTU_db.fasta",
                     "--all_fasta", "dada2_seq_and_OTU_db.fasta",
                     "--seq_db", OTU_DATABASE,
                     "--min_novel_cluster_threshold",
                     "0",
                     "--left", LEFT_READS,
                     "--right", RIGHT_READS,
                     "--Name_of_project",
                     os.path.join(DADA2_FOLDER, "clusters"),
                     "--in",
                     SWARM_OUT,
                     "--difference", str(SWARM_D_VALUE),
                     "-o",
                     os.path.join(DADA2_FOLDER,
                                  "%s_swarm_results_%d.RESULTS" %
                                  (PREFIX, SWARM_D_VALUE))]
            swarm_result = os.path.join(DADA2_FOLDER,
                                        "%s_swarm_results_%d.RESULTS" %
                                        (PREFIX, SWARM_D_VALUE))
            RESULTS.append("swarm\t%s" % swarm_result)
            cmd_s = ' '.join(cmd_s)
            if args.align:
                logger.info("going to align the cluster. Will take ages!")
                cmd_s = cmd_s + " --align True"
            if args.percent_identity:
                cmd_s = cmd_s + " --blast True"
            logger.info("%s = post analysis command", cmd_s)
            pipe = subprocess.run(cmd_s, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True)
            logger.info("graphically represent swarm clusters")
            plot_cmd = ["python",
                        os.path.join(FILE_DIRECTORY,
                                     "bin",
                                     "draw_bar_chart_of_clusters.py"),
                        "-i",
                        SWARM_OUT,
                        " --db",
                        OTU_DATABASE_SWARM]
            plot_cmd = ' '.join(plot_cmd)
            logger.info("plotting command = %s", plot_cmd)
            pipe = subprocess.run(plot_cmd, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True)
            # move the dada2 fasta and Rplots to the folder
            subprocess.run(["mv", PREFIX + "_DADA2.fasta",
                            os.path.join(PREFIX + "_RESULTS",
                                         PREFIX + "_DADA2.fasta")])
            subprocess.run(["mv", "Rplots.pdf",
                            os.path.join(PREFIX + "_RESULTS",
                                         "Rplots.pdf")])

    #####################################################################
    # run BLASTCLUST
##    if "blastclust" in tools_list:
##        print("we are not going to run blastclust")
##        continue  # this is just too slow
##        Run_blastclust = args.Run_blastclust.upper()
##        if Run_blastclust.rstrip() == "YES":
##            blastclust_threshold = args.blastclust_threshold  # for now
##            bc = PREFIX + "_blastclust_%s" % str(blastclust_threshold)
##            BLASTCL_FOLDER = make_folder(os.path.join(RESULT_FOLDER,
##                                                      bc),
##                                         WORKING_DIR)
##            logger.info("running blastclust. This is slow")
##            bc = blast.Blastclust("blastclust")
##            result = bc.run("assembled_fa_and_OTU_db.fasta", BLASTCL_FOLDER,
##                            THREADS)
##
##            # reformat_swarm_cls(swarm, db, db_and_reads, outfile)
##            # use this to reformat the blastclust clusters
##            result_file_r = (os.path.join(BLASTCL_FOLDER,
##                                          "assembled_fa_and_OTU_db.fasta" +
##                                          ".blastclust99.forR"))
##            reformat_swarm_cls(os.path.join(BLASTCL_FOLDER,
##                                            "assembled_fa_and_OTU_db.fasta" +
##                                            ".blastclust99.lst"),
##                               OTU_DATABASE,
##                               "assembled_fa_and_OTU_db.fasta",
##                               result_file_r,
##                               False)
##            # add this file for Rand index comparison later
##            CLUSTER_FILES_FOR_RAND_INDEX.append(result_file_r)
##            logger.info("graphically represent swarm clusters")
##            plot_cmd = ["python",
##                        os.path.join(FILE_DIRECTORY,
##                                     "bin",
##                                     "draw_bar_chart_of_clusters.py"),
##                        "-i",
##                        result.outfilename,
##                        " --db",
##                        OTU_DATABASE]
##            plot_cmd = ' '.join(plot_cmd)
##            logger.info("plotting command = %s", plot_cmd)
##            pipe = subprocess.run(plot_cmd, shell=True,
##                                  stdout=subprocess.PIPE,
##                                  stderr=subprocess.PIPE,
##                                  check=True)
##            cmd_b = ["python",
##                     os.path.join(FILE_DIRECTORY,
##                                  "post_analysis",
##                                  "get_results_from_cluster_and_novel_" +
##                                  "clusterings_cd_hit.py"),
##                     "-f", ASSEMBLED + ".bio.chopped.fasta",
##                     "--all_fasta", "assembled_fa_and_OTU_db.fasta",
##                     "--seq_db", OTU_DATABASE,
##                     "--min_novel_cluster_threshold",
##                     args.min_novel_cluster_threshold,
##                     "--left", LEFT_READS,
##                     "--right", RIGHT_READS,
##                     "--Name_of_project",
##                     os.path.join(BLASTCL_FOLDER, "clusters_%s" %
##                                  str(blastclust_threshold)),
##                     "--in",
##                     result.outfilename,
##                     "--difference", str(blastclust_threshold),
##                     "-o",
##                     os.path.join(BLASTCL_FOLDER,
##                                  "%s_blastclust_results_%s.RESULTS" %
##                                  (PREFIX, str(blastclust_threshold))),
##                     "--old_to_new", "db_old_to_new_names.txt"]
##
##            bc_result = os.path.join(BLASTCL_FOLDER,
##                                     "%s_blastclust_results_%s.RESULTS" %
##                                     (PREFIX, str(blastclust_threshold)))
##            RESULTS.append("blastclust\t%s" % bc_result)
##
##            cmd_b = ' '.join(cmd_b)
##            if args.align:
##                logger.info("going to align the cluster. Will take ages!")
##                cmd_b = cmd_b + " --align True"
##            if args.percent_identity:
##                cmd_b = cmd_b + " --blast True"
##            logger.info("%s = post analysis command", cmd_b)
##            pipe = subprocess.run(cmd_b, shell=True,
##                                  stdout=subprocess.PIPE,
##                                  stderr=subprocess.PIPE,
##                                  check=True)
    try:
        Rand_results = pairwise_comparison_Rand(CLUSTER_FILES_FOR_RAND_INDEX,
                                                os.path.join(RESULT_FOLDER,
                                                             PREFIX +
                                                             "_Rand_compar" +
                                                             "ison.txt"))
        #for comp in Rand_results:
            #logger.info("Rand comparison: %s", comp)
    except ValueError:
        logger.warning("Rand comparison failed.")
    # compress the reads to save space
    compress(LEFT_READS)
    compress(RIGHT_READS)
    logger.info("compressed the reads")

    #####################################################################
    # compare the results files
    compare_prog = os.path.join(FILE_DIRECTORY,
                                "bin",
                                "compare_results.py")
    blastclust_threshold = str(args.blastclust_threshold)
    comp = "%s_RESULTS_cd_%s_sw_%s_BC_%s_V_%s.txt" % (PREFIX,
                                                      CDHIT_THRESHOLD,
                                                      str(SWARM_D_VALUE),
                                                      blastclust_threshold,
                                                      str(VSEARCH_THRESHOLD))
    # write the result file name to file. Easier to get these for the
    # next part - compare_prog
    f_out = open("temp.txt", "w")
    for result in RESULTS:
        f_out.write(result + "\n")
    f_out.close()
    cmd_r = " ".join(["python",
                      compare_prog,
                      " -o",
                      os.path.join(RESULT_FOLDER, comp),
                      " --in_list",
                      "temp.txt"])
    logger.info("%s = comparison comment", cmd_r)
    try:
        pipe = subprocess.run(cmd_r, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
    except ValueError:
        logger.warning("overall comparison failed.")
    if args.cleanup.upper() == "YES":
        # remove loads of files for the user
        # as previously having other files confused people
        remove_list = [ASSEMBLED + "for_swarm.fasta",
                       ASSEMBLED + ".bio.chopped.fasta",
                       OTU_DATABASE_SWARM,
                       SWARM_OUT,
                       "temp.fasta",
                       "temp.txt",
                       "assembled_fa_and_OTU_db.fasta",
                       "assembled_reads_and_OTU_db.fasta",
                       "db_old_to_new_names.txt",
                       "db_old_to_new_names_vsearch.txt",
                       "assembled_fa_and_OTU_db_vesearch.fasta",
                       db_derep,
                       (os.path.join(TRIM_FOLDER,
                                     PREFIX + "_unpaired_R1.fq.gz")),
                       (os.path.join(TRIM_FOLDER,
                                     PREFIX + "_unpaired_R2.fq.gz")),
                       ASSEMBLED + "drep.vsearch.fasta",
                       "OTU.1.bt2",
                       "OTU.2.bt2",
                       "OTU.3.bt2",
                       "OTU.4.bt2",
                       "OTU.rev.1.bt2",
                       "OTU.rev.2.bt2",
                       "error.log",
                       "dada2_seq_and_OTU_db.fasta",
                       "dada2.R",
                       results_pear.outfilediscarded,
                       results_pear.outfileunassmbledfwd,
                       results_pear.outfileunassembledrev,
                       PREFIX + "_R1_001_primers_trimmed.fastq.gz",
                       PREFIX + "_R2_001_primers_trimmed.fastq.gz",
                       "sequence_table.txt",
                       "temp",
                       results_pear.outfilediscarded,
                       results_pear.outfileunassmbledfwd,
                       results_pear.outfileunassembledrev,
                       results_pear.outfileassembled]
        for unwanted in remove_list:
            try:
                os.remove(unwanted)
                logger.info("deleting: %s", unwanted)
            except:
                logger.info("could not find %s", unwanted)
    shutil.rmtree(PREFIX)
    if args.Run_dada2.upper() == "YES":
        shutil.rmtree("dada2")
        shutil.rmtree("filtered")
    if ERROR_CORRECTION and args.cleanup:
        # this folders will only be there is EC was run
        shutil.rmtree(os.path.join(EC_FOLDER))
    logger.info("Pipeline complete: %s", time.asctime())
