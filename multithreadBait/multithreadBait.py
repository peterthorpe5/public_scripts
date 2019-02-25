#!/usr/bin/env python3
# title: run multi thread instances of mirabait
# Author: Peter Thorpe. 
# why? speed"""

from threading import Thread
import time
import logging
import logging.handlers
import subprocess
import sys
import os
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator


usage = """Use as follows:

script to run mirabait on multi threads.

This will only split at num of scaffold. It will NOT split up scaffold as
these could interefer with real hits.

requires: Biopython
"""

# check the user is using python 3. 
if sys.version_info[:1] != (3,):
    # e.g. sys.version_info(major=3, minor=5, micro=2,
    # releaselevel='final', serial=0)
    # break the program
    print ("currently using:", sys.version_info,
           "  version of python")
    raise ImportError("Python 3.x is required for TranStart.py")
    print ("if you want to force it to use python 2 replace this line"
           " sys.version_info[:1] != (3) with "
           " sys.version_info[:1] != (2) ")
    sys.exit(1)


######
def get_args():
    parser = argparse.ArgumentParser(description="multi core " +
                                     "mirabait ",
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("multi")[0]
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-k", "--kmer", dest="kmer",
                      default="65",
                      help="kmer length for mira bait match")

    optional.add_argument("-o", "--out",
                          dest="outfile",
                          default=None,
                          help="Output filename",
                          metavar="FILE")

    optional.add_argument("-g", "--genome",
                          dest="genome",
                          default=None,
                          help="the genome sequence",
                          metavar="FILE")

    optional.add_argument("-b", "--bait",
                          dest="bait",
                          default=None,
                          help="bait to do the kmer matching with",
                          metavar="FILE")

    optional.add_argument("-l", "--left",
                          dest="left",
                          default=None,
                          help="left reads",
                          metavar="FILE")

    optional.add_argument("-r", "--right",
                          dest="right",
                          default=None,
                          help="right reads",
                          metavar="FILE")

    optional.add_argument("-t", "--threads",
                          dest="threads",
                          default=None,
                          help="number of threads to use")

    optional.add_argument("--remove_file",
                          dest="remove_file",
                          default="YES",
                          help="clean up of not: yes or no")


    optional.add_argument("--logfile",
                          dest="logfile",
                          action="store",
                          default="pipeline.log",
                          type=str,
                          help="Logfile name")
    args = parser.parse_args()
    return args, file_directory


def subproces_func(cmd):
    """func to run subporcess"""
    print("running ...", cmd)
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
    return pipe


def make_me_a_folder(infolder):
    directory = os.getcwd()
    try:
        os.mkdir(infolder)
    except OSError:
        print("already exists")
    #os.chdir("split_fasta_files")


def decompress(infile, threads):
    """function to decompress gzipped reads"""
    cmd = ' '.join(["pigz", "-d", "--quiet", "--force",
                    "-p", str(threads), infile])
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
    return os.path.splitext(infile)[0]


def compress(infile, threads):
    """function to compress reads"""
    cmd = ' '.join(["pigz", "--quiet", "--force",
                    "-p", str(threads), infile])
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)


def test_reads_exist_and_suffix(reads):
    """function to test if the reads exist and if there
    are compressed or not.
    Take in a read file.
    This is required, if a run is cancelled then re-run
    with the same command and the reads have been decompresssed,
    it will fail.
    return the proper name the reads will have"""
    if not os.path.isfile(reads):
        # there is something wrong ...
        # this file doesnt exist. Maybe already decompressed
        # but not what the user is asking for
        if os.path.isfile(os.path.splitext(reads)[0]):
            # the reads have already been decompressed
            return os.path.splitext(reads)[0]
    if os.path.isfile(reads):  # file is real
        if reads.endswith(".gz"):  # decompress them
            return reads
    else:
        error = "\nERROR: %s   FILE DOES NOT EXIT" % reads
        sys.exit(error + ". Check your input file path and name\n")


def count_fq_reads(in_fastq):
    """function count the number of reads"""
    # open the fastq file
    in_file = open(in_fastq)
    # iterate through the fastq file
    total_reads = 0
    for (title, sequence, quality) in \
        FastqGeneralIterator(in_file):
        total_reads = total_reads+1
    in_file.close()
    return total_reads


def out_fq_name(out_fastq, count):
    outname = out_fastq.split(".f")[0] + "_" + str(count) + ".fastq"
    outfile = os.path.join("temp_reads", outname)
    return outfile


def write_out_fq (in_fastq, threads, number_of_seq, out_fastq):
    #open the fastq file
    number_of_seq = int(number_of_seq)
    in_file = open(in_fastq)
    # enumerate is a way of counting i
    # iterate through the fastq file
    list_of_fq_file = set([])
    count = 1
    outfq = out_fq_name(out_fastq, count)
    out_file = open(outfq, "w")
    list_of_fq_file.add(outfq)

    for i, (title, seq, qual) in enumerate(FastqGeneralIterator(in_file)):
        #python magic to identify every number_of_seq "loop"
        out_file.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

        if i % number_of_seq == 0:
            #write this to a file
            out_file.close()
            count = count + 1
            outfq = out_fq_name(out_fastq, count)
            out_file = open(outfq, "w")
            list_of_fq_file.add(outfq)            
    in_file.close()
    return list_of_fq_file



#######################################################################
# Run as script
if __name__ == '__main__':
    args, FILE_DIRECTORY = get_args()
    # Set up logging
    logger = logging.getLogger('multithreadbait.py: %s' % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)
    if args.logfile == "pipeline.log":
        args.logfile = "mirabait_multithread" + "_pipeline.log"
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
    PREFIX = args.left.split("R")[0]
    logger.info(sys.version_info)
    logger.info("Command-line: %s", ' '.join(sys.argv))
    logger.info("Starting testing: %s", time.asctime())
    THREADS = args.threads
    make_me_a_folder("temp_reads")
    # test to see if they have been decompressed already
    args.left = test_reads_exist_and_suffix(args.left)
    args.right = test_reads_exist_and_suffix(args.right)
    if args.left.endswith(".gz"):
        logger.info("decompressing left reads")
        args.left = decompress(args.left, args.threads)
        logger.info(args.left)
    if args.right.endswith(".gz"):
        logger.info("decompressing right reads")
        args.right = decompress(args.right, args.threads)
    logger.info("counting number of reads")
    total_reads = count_fq_reads(args.left)
    logger.info("number of reads = %d", total_reads)
    # number of file will be determined by number of threads - 1
    NumSeqPerOutfiles = (int(total_reads)/int(THREADS))-1
    logger.info("number of reads per file = %d", NumSeqPerOutfiles)
    # write_out_fq(in_fastq, threads, number_of_seq, out_fastq)
    left_fq_file = write_out_fq(args.left, THREADS,
                                NumSeqPerOutfiles,
                                args.left)
    right_fq_file = write_out_fq(args.right, THREADS,
                                 NumSeqPerOutfiles,
                                 args.right)
    # start of multi threading
    # #/mirabait -I -k 99 -b clc.all.fa.contim_filtered.fasta
    # -p R1.fastq R2.fastq
    mira_cmds = set([])
    for reads in left_fq_file:
        left_test_reads = reads
        right_test_reads = reads.replace("R1", "R2")
        cmd = " ".join(["mirabait",
                        "-I",
                        "-k %s" % args.kmer,
                        "-b %s" % args.bait,
                        "-p",
                        left_test_reads,
                        right_test_reads])
        mira_cmds.add(cmd)
    # split this over all the cores asked for
    threads = []
    for i in mira_cmds:
        logger.info(i)
        t = Thread(target=subproces_func, args=(i,))
        threads.append(t)
        t.start()

    # Wait for all of them to finish
    for x in threads:
        x.join()
        
    logger.info("compressing original input files")
    compress(args.left, THREADS)
    compress(args.right, THREADS)
    
    make_me_a_folder("match")
    make_me_a_folder("miss")
    #
    outfile = os.path.join("match" , "bait_match_" + PREFIX + "R1.fastq")
    cmd = " ".join(["cat", "bait_match_*R1*",
                    ">",
                    outfile])
    logger.info(cmd)

    subproces_func(cmd)
    logger.info("compressing baited match R1 reads")
    compress(outfile, THREADS)
    print("compressed")
    #
    outfile = os.path.join("miss" , "bait_miss_" + PREFIX + "R1.fastq")
    cmd = " ".join(["cat", "bait_miss_*R1*",
                    ">",
                    outfile])
    logger.info(cmd)

    subproces_func(cmd)
    compress(outfile, THREADS)

    #
    outfile = os.path.join("match" , "bait_match_" + PREFIX + "R2.fastq")
    cmd = " ".join(["cat", "bait_match_*R2*",
                    ">",
                    outfile])
    subproces_func(cmd)
    compress(outfile, THREADS)

    outfile = os.path.join("miss" , "bait_miss_" + PREFIX + "R2.fastq")
    cmd = " ".join(["cat", "bait_miss_*R2*",
                    ">",
                    outfile])
    subproces_func(cmd)
    compress(outfile, THREADS)

    if "Y" in args.remove_file.upper():
        logger.info("removing files")
        cmd = " ".join(["rm ", "bait_m*_*"])
        subproces_func(cmd)
        cmd = " ".join(["rm", "-rf", "temp_reads"])
        subproces_func(cmd)



    subproces_func(cmd)
    



        
