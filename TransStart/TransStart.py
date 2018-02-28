#!/usr/bin/env python
# title: Find the start of transcription based on per base read count
# (c) The James Hutton Institute 2018
# Author: Peter Thorpe. The James Hutton Institute, Uk
# why? Start of transcription is interesting! THIS IS NOT finding the start
# of coding squence. but where transcription binds and starts

# imports
import subprocess
import sys
import os
from optparse import OptionParser
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy
import datetime
import logging
import logging.handlers
import time
# TODO make re look for Kozak GCC[A/G] CCatg
import re


if "-v" in sys.argv or "--version" in sys.argv:
    print("Transcription start finder based on read depth coverage v0.1.0")
    cmd = "samtools 2>&1 | grep -i ^Version"
    sys.exit(os.system(cmd))

if "--help_full" or "-h" in sys.argv:
    print("""Transcription start finder based on read depth coverage v0.1.0
 why? Start of transcription is interesting! THIS IS NOT finding the start
# of coding squence. but where transcription binds and starts

Requires:
samtools 1.2 or greater
Biopython
numpy

steps:

1) Index your transcriptome and map your reads back to your genome:
    Use STAR: a spice aware aligner
2) STAR --runMode genomeGenerate --runThreadN 10 --limitGenomeGenerateRAM
    74554136874 --genomeDir /PATH_TO/M.cerasi/star_indicies
    --genomeFastaFiles /PATH_TO/M.cerasi/Mc_v1.fasta
3) STAR --genomeDir star_indicies/ --runThreadN 12
    --outFilterMultimapNmax 5 --outSAMtype BAM SortedByCoordinate
    --outFilterMismatchNmax 7 --readFilesCommand zcat
    --outFileNamePrefix Mc --readFilesIn /PATH_TO/M.cerasi/RNAseq/Mc_R1.fq.gz
    /PATH_TO/M.cerasi/RNAseq/Mc_R2.fq.gz
2) Index your genome:
4) sort your bam out!
    samtools -@ 12 sort unsorted.bam sorted.bam
5) Index your sorted.bam
    samtools index sorted.bam

How?

This script looks at the number of reads that map per base for your gene
of interest.
If the number of mapped reads is greater than "threshold (default=3, reads
mapped per current base)" then it performs statistical
analysis on this to determine if the start of transcription
    """)
    cmd = "samtools 2>&1 | grep -i ^Version"
    sys.exit(os.system(cmd))



# make a temp_folder_for_all_the_out_files
dest_dir = os.path.join(os.getcwd(),
                        'temp_fix_five_prime')
try:
    os.makedirs(dest_dir)
except OSError:
    print ("folder already exists, I will write over what is in there!!")

##############################################################################
# functions

def sys.exit(msg, error_level=1):
   """Print error message to stdout and quit with given error level."""
   sys.stderr.write("%s\n" % msg)
   sys.exit(error_level)


def get_total_coverage(bam_file, outfile):
    """ function to get the total coverage a found in the bam file
    Takes in a bam file. Call samtools idxstats to obtain values
    Results are put in ./temp_fix_five_prime/ and written over
    each time
    funtion returns a dictions off overall expression
    returns a dictionary: key[transcript], vals = ['577', '274', '0'] len,
    # reads_mapped, last_coloumn
    """
    # Run samtools idxstats (this get the coverage for all transcripts:
    # assigne the outfile with the temp folder to keep thing more tidy
    oufile_dir_file = os.path.join("temp_fix_five_prime",
                                   outfile)
    cmd = " ".join(['samtools',
                    'idxstats',
                    bam_file,
                    '>',
                    oufile_dir_file])
    # data was saved in idxstats_filename
    # call the func
    pipe = subproces_func(cmd)
    # creat a dictioanry to hold all the total expression values for the
    # transcripts.
    overall_expression_dic = dict()
    with open(oufile_dir_file, "r") as handle:
        for line in handle:
            data = line.rstrip("\n").split("\t")
            transcript = data[0]
            overall_expression_dic[transcript] = [int(x) for x in data[1:]]
    # print overall_expression_dic["Mp_O_20647_c0_seq2"]
    # returns a dictionary: key[transcript], vals = ['577', '274', '0'] len,
    # reads_mapped, last_coloumn
    return overall_expression_dic


def index_genome_file(genome):
    """index the genome file
    return a seq_record object that can be accessed
    in a dictionary like manner"""
    genome_database = SeqIO.index(genome, "fasta")
    return genome_database


def index_gff(gff):
    """take in gff file. Only looks at gene features. Returns many dicts"""
    f_in = open(gff, "r")
    gene_start_stop_dict = dict()
    gene_scaff_dict = dict()
    gene_first_exon_dict = dict()
    gene_direction = dict()
    gene_gff_line = dict()
    gene_set = set([])
    for line in f_in:
        if line.startswith("#"):
            continue
        if not line.strip():
            continue
        assert len(line.split("\t")) == 9 ,"GFF fields wrong length should be 9"
        scaff, source, feature, start, stop, score, \
              direction, frame, gene_info = line.split("\t")
        gene_info = gene_info.replace("ID=", "").split()[0]
        gene_info = gene_info.split(".t")[0]
        gene_info = gene_info.replace(";", "")
        gene_info = gene_info.split("Note=")[0]
        if not gene_first_exon_dict[gene]:
            if feature == "exon":
                start_stop = "%s\t%s" % (start, stop)
                gene_first_exon_dict[gene] = start_stop
        if not feature == "gene":
            continue
        gene_gff_line[gene] = line
        gene_set.add(gene_info)
        start_stop = "%s\t%s" % (start, stop)
        gene_start_stop_dict[gene_info] = start_stop
        gene_scaff_dict[gene] = scaff
        gene_direction[gene] = direction
    f_in.close()
    return gene_start_stop_dict, gene_first_exon_dict, \
           gene_scaff_dict, gene_direction, gene_set, gene_gff_line


def average_standard_dev(positions):
    """function to return the avaerage and stadard deviation
    for a list of number.
    Uses Numpy to do the calculations"""
    the_mean = numpy.mean(positions)
    standard_dev = numpy.std(positions)
    return the_mean, standard_dev


def mean_coverage(coverage_array, slice_start, slice_end):
    """function to get the mean coverage for a sliced region"""
    selected_coverage = coverage_array[slice_start : slice_end]
    return mean(selected_coverage)


def middle_portion_of_transcript(seq):
    """return coordinates for middle portion"""
    coord_lower = int(len(seq)/4.0)
    coord_upper = int(len(seq)-coord_lower)
    return coord_lower, coord_upper


def subproces_func(cmd):
    """func to run subporcess"""
    pipe = (cmd, shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True)
    return pipe


def fill_in_zero_cov(all_coverage, depth_filename):
    """ function to fil in the zero values returned by samtools depth,
    or not returned by"""
    for line in open(depth_filename):
        ref, possition, coverage = line.rstrip("\n").split("\t")
        possition = int(possition) - 1
        # assign the correct coverage to the relavent postiton,
        # thus removing the zero value.
        # Or if there is no info, it remains as a zero.
        all_coverage[possition] = int(coverage)
    return depth_filename, all_coverage


def run_samtools_depth(scaffold_start_stop, bam_file, outfile):
    """funct to run samtools depth"""
    cmd = " ".join(["samtools",
                    "depth",
                    "-r",
                    scaffold_start_stop,
                    bam_file,
                    ">",
                    outfile])
    # call the func to run subprocess
    pipe = subproces_func(cmd)
    return pipe


def walk_away_from_start(start, stop,
                         direction, walk=3):
    """function to walk away from the start to get the coverage stats"""
    if direction == "+":
        current_start = int(start) - int(walk)
        current_end = int(start)
    else:
        assert direction == "-", "direction does sign error!"
        current_end = stop + int(walk)
        current_start = int(stop)
    return current_start, current_stop


def add_one_on_direction_aware(current_start, current_end, interation_value,
                               direction):
    """function to alter the current start and end values based on the direction"""
    if direction == "+":
        current_start = current_start - interation_value
        current_end = current_end - interation_value
    else:
        current_start = current_start + interation_value
        current_end = current_end + interation_value
    return current_start, current_end


##### main function
def TranscriptionFind(genome,
                      gene_start_stop_dict,
                      gene_first_exon_dict,
                      gene_scaff_dict,
                      gene_direction,
                      gene_set,
                      gene_gff_line,
                      bam,
                      stand_dev_threshold,
                      walk,
                      min_value,
                      interation_value,
                      out_file):
    """iterate through gene in gff. Call samtools, get expression
    Does the expression fall within X standard deviations when compare
    to the first exon?
    """
    # open outfile:
    interation_value = int(interation_value)
    file_out = open(out_file, "w")
    out_str = "\t".join(["#gene",
                         "transcription start based on stats:",
                         "start",
                         "stop",
                         "transcription start based on min value:",
                         "start",
                         "stop",
                         "gene_start",
                         "gene_stop"
                         "gene_mean_RNAseq_cov",
                         "gene_std_RNAseq_cov",
                         "for 1st exon",
                         "exon mean",
                         "exom std",
                         "coding direction\n"])
    file_out.write(out_str)

    # threshold for min_read_count
    min_read_count_threshold = int(min_read_count) # default is 30

    for gene in gene_set:
        gene = gene.rstrip()
        start, stop = gene_start_stop_dict[gene]
        scaffold = gene_scaff_dict[gene]
        direction = gene_direction[gene]
        exon_start, exon_stop = gene_first_exon_dict[gene]

        # call samtools to get the depth per posititon for
        # the transcript of interest
        depth_filename = os.path.join("temp_fix_five_prime",
                                      "depth.tmp")
        exon_depth_file = os.path.join("temp_fix_five_prime",
                                      "exon_depth.tmp")
        scaffold_start_stop = "%s:%s-%s" %(scaffold, start, stop)
        # call the func to run
        pipe = run_samtools_depth(scaffold_start_stop, bam_file,
                                  depth_filename)
        scaffold_exon_start_stop = "%s:%s-%s" %(scaffold, exon_start,
                                                exon_stop)
        pipe = run_samtools_depth(scaffold_exon_start_stop, bam_file,
                                  exon_depth_file)
        # assign zeros to all positions of the transcript,
        # as samtool does no report zeros
        all_coverage = [0] * len(scaffold)
        all_coverage, depth_filename = fill_in_zero_cov(all_coverage, depth_filename)
        exon_all_coverage, exon_depth_file = fill_in_zero_cov(all_coverage,
                                                              exon_depth_file)
        # get the mean and std reads per base for exon 1
        exon_mean, exon_stdDev = average_standard_dev(exon_all_coverage
                                                     [exon_start:exon_stop])
        # get the mean and std reads per base for exon 1
        gene_mean, gene_stdDev = average_standard_dev(all_coverage
                                                     [start:stop])
        out_str = "\t".join([gene + ":",
                            "Cov min: %i" % min(all_coverage),
                            "max: %i" % max(all_coverage),
                            "gene mean %0.2f:" % gene_mean,
                            "gene std %0.2f:" % gene_stdDev,
                            "Sliced section:",
                            "exon mean %0.2f" % exon_mean,
                            "exom std: %0.2f" % exon_stdDev,
                            direction,
                            "\n"])
        logger.info(out_str)

        cut_off = exon_mean - (int(stand_dev_threshold)
                               * exon_stdDev)

        position_mean_cov = mean(all_coverage[exon_start:exon_stop])
        # walk in 3 bases to find the position where coverage sig drops
        current_end = stop
        current_start = start
        while position_mean_cov < cut_off:
            current_start, current_end = walk_away_from_start(current_start,
                                                              current_end,
                                                              direction, walk)
            current_start, current_end = add_one_on_direction_aware(current_start,
                                                                    current_end,
                                                                    interation_value,
                                                                    direction)
            position_mean_cov = mean(all_coverage[current_start:current_end])

        # run the while loop again to find the position where the expression
        # is less than the option min value
        current_end1 = stop
        current_start1 = start
        while position_mean_cov < int(min_value):
            current_start1, current_end1 = walk_away_from_start(current_start1,
                                                                current_end1,
                                                                direction, walk)
            current_start1, current_end1 = add_one_on_direction_aware(current_start1,
                                                                      current_end1,
                                                                      interation_value,
                                                                      direction)
            position_mean_cov = mean(all_coverage[current_start1:current_end1])

        out_str = "\t".join([gene + ":",
                             "transcription start based on stats:",
                             str(current_start),
                             str(current_stop),
                             "transcription start based on min value:",
                             str(current_start1),
                             str(current_stop1),

                             "gene start %s" % start,
                             "gene ened %s" % stop,
                             "gene mean %0.2f:" % gene_mean,
                             "gene std %0.2f:" % gene_stdDev,
                             "Sliced section:",
                             "exon mean %0.2f" % exon_mean,
                             "exom std: %0.2f" % exon_stdDev,
                             direction,
                             "\n"])
        logger.info("main function finished")


###############################################################################################
#to run it:


usage = """Use as follows:

Transcription start finder based on read depth coverage v0.1.0
THIS IS NOT finding the start
# of coding squence. but where transcription binds and starts

python TranscriptionStart.py -g genome.fasta
    --bam index_sorted_bam_file.bam
    --gff genes.gff -o outfile


Requirements:
    you must have samtools in your PATH
    Biopython
    numpy


steps:

1) Map RNAseq using STAR
2) sort your bam out!
    samtools sort unsorted.bam sorted.bam
3) Index your sorted.bam
    samtools index sorted.bam

"""



parser = OptionParser(usage=usage)


parser.add_option("--gff", dest="gff",
                  default=None,
                  help="the predicted coordinate for the cds predictions .gff",
                  metavar="FILE")

parser.add_option("-g", "--genome", dest="genome",
                  default=None,
                  help="the genome sequence. Not currently used. TO DO",
                  metavar="FILE")

parser.add_option("--bam", dest="bam",
                  default=None,
                  help="the sorted, indexed bam file as a result of " +
                  "the reads being mapped back to the transcriptome  .bam",
                  metavar="FILE")

parser.add_option("--std",
                  dest="stand_dev_threshold",
                  default=5,
                  help="If the expression of the current region is less then "
                  " this number of std away from the first exon mean. Default = 5")

parser.add_option("--walk",
                  dest="walk",
                  default=3,
                  help="the number of bases to walk away to find expression " +
                  "Default = 3")

parser.add_option("--min_value",
                  dest="min_value",
                  default=2,
                  help="the min expression to return the coordinates for " +
                  "Default = 2")

parser.add_option("--interation_value",
                  dest="interation_value",
                  default=1,
                  help="size of window to walk " +
                  "Default = 1")

parser.add_option("--help_full", dest="help_full", default=None,
                  help="prints out a full description of this program")

parser.add_option("--logger", dest="logger", default=False,
                  help="Output logger filename. Default: " +
                  "outfile_std.log",
                  metavar="FILE")

parser.add_option("-o", "--out", dest="outfile", default="results.out",
                  help="Output filename (default: results.out)",
                  metavar="FILE")

(options, args) = parser.parse_args()


#-g
genome = options.genome
#--gff
gff = options.gff
#--bam
bam_file = options.bam
#-o
outfile= options.outfile
# ==stand_dev_threshold
stand_dev_threshold = options.stand_dev_threshold
# logfile
logger = options.logger
if not logger:
    log_out = "%s_%s_std.log" % (outfile, stand_dev_threshold)

if not os.path.isfile(transcriptome_file):
    sys.exit("Input transcriptome file not found: %s" % transcriptome)

if not os.path.isfile(bam_file):
    sys.exit("Input BAM file not found: %s" % bam)

if not os.path.isfile(bam_file + ".bai"):
    print("you have not indexed you bam file. I will do it.")
    cmd = " ".join(["samtools",
                    "index",
                    bam_file])
    # call the func
    pipe = subproces_func(cmd)


#######################################################################
# Run as script
if __name__ == '__main__':
    if not os.path.isfile(genome):
        sys.exit("Input genome file not found: %s" % genome)
    if not os.path.isfile(genes_fasta):
        sys.exit("Input genes_fasta file not found: %s" % genes_fasta)
    if not os.path.isfile(bam_file):
        sys.exit("Input BAM file not found: %s" % bam)

    # Set up logging
    logger = logging.getLogger('Transcription_start_predictor.py: %s' % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)
    try:
        logging.basicConfig(filename='example.log',level=logging.DEBUG)
        logstream = open(log_out, 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        logger.error("Could not open %s for logging", logger)
        sys.exit(1)
    # Report input arguments
    logger.info(sys.version_info)
    logger.info("Command-line: %s", ' '.join(sys.argv))
    logger.info("Starting testing: %s", time.asctime())
    logger.info("get total expression")
    overall_expression_dic = get_total_coverage(bam_file,
                                                over_all_expression)
    logger.info("Indexidng gff: ", gff)
    gene_start_stop_dict, gene_first_exon_dict,\
                          gene_scaff_dict, gene_direction, gene_set,\
                          gene_gff_line = index_gff(gff)
    logger.info("runnin analysis")
    TranscriptionFind(genome,
                      gene_start_stop_dict,
                      gene_first_exon_dict,
                      gene_scaff_dict,
                      gene_direction,
                      gene_set,
                      gene_gff_line,
                      bam,
                      stand_dev_threshold,
                      option.walk,
                      option.min_value,
                      option.interation_value,
                      out_file)
    logger.info("Pipeline complete: %s", time.asctime())

