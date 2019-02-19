#!/usr/bin/env python3
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

if "-v" in sys.argv or "--version" in sys.argv:
    print("TransStart.py Transcription start finder "
          "based on read depth coverage v0.1.0")
    cmd = "samtools --version"
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
    sys.exit("samtool version installed = ", pipe)

if "--help_full" in sys.argv:
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
    --outFileNamePrefix Mc --readFilesIn /RNAseq/Mc_R1.fq.gz
    /RNAseq/Mc_R2.fq.gz
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


 python TransStart.py -g genome.fasta
    --bam index_sorted_bam_file.bam
    --walk 5 --interation_value 1 
    --gff genes.gff -o outfile


Requirements:
    you must have in your PATH:
    samtools
    Biopython
    numpy

    """)
    cmd = "samtools --version"
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
    sys.exit("samtool version installed = ", pipe)

##############################################################################
# functions

def get_total_coverage(bam_file, outfile):
    """ function to get the total coverage a found in the bam file
    Takes in a bam file. Call samtools idxstats to obtain values
    Results are put in ./temp_reads_per_base/ and written over
    each time
    funtion returns a dictions off overall expression
    returns a dictionary: key[transcript], vals = ['577', '274', '0'] len,
    # reads_mapped, last_coloumn
    """
    # Run samtools idxstats (this get the coverage for all transcripts:
    # assigne the outfile with the temp folder to keep thing more tidy
    oufile_dir_file = os.path.join("temp_reads_per_base",
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


def index_genome_file(genome, logger):
    """index the genome file
    return a seq_record object that can be accessed
    in a dictionary like manner"""
    logger.info("indexing genome: %s", genome)
    genome_database = SeqIO.index(genome, "fasta")
    logger.info("indexing genome finished")
    return genome_database


def split_gene_name(gene):
    """dirty funk to clean up the gene name.
    Just want to say. I hate all the different ways gene names
    can be messed around with in gff files"""
    gene_info = gene.replace("ID=", "").split()[0]
    gene_info = gene_info.split(".t")[0]
    if "-T" in gene_info:
            gene_info = gene_info.split("-T")[0] # funannotate models
    gene_info = gene_info.replace(";", "")
    gene_info = gene_info.replace("Parent=", "")
    gene_info = gene_info.split("Note=")[0]
    gene_info = gene_info.split("Name=")[0]
    return gene_info


def index_gff(gff, logger):
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
        assert len(line.split("\t")) == 9 , "GFF fields wrong length should be 9"
        scaff, source, feature, start, stop, score, \
              direction, frame, gene_info = line.split("\t")
        gene = split_gene_name(gene_info)
        scaff = scaff.rstrip()
        if feature == "gene":
            gene_gff_line[gene] = line
            gene_set.add(gene)
            start_stop = "%s\t%s" % (start, stop)
            gene_start_stop_dict[gene] = start_stop
            gene_scaff_dict[gene] = scaff
            gene_direction[gene] = direction
        if not gene in gene_first_exon_dict.keys():
            if feature == "exon" or feature == "CDS":
                start_stop = "%s\t%s" % (start, stop)
                gene_first_exon_dict[gene] = start_stop
    f_in.close()
    logger.info("Number of genes = %d", len(gene_set))
    return gene_start_stop_dict, gene_first_exon_dict, \
           gene_scaff_dict, gene_direction, gene_set, gene_gff_line


def create_gff_line(gffline, gene, TranStart, TranStop):
    """func to create a gff line with a new name and
    the start of transcription and end of transcription"""
    scaff, source, feature, start, stop, score, \
              direction, frame, gene_info = gffline.split("\t")
    gene_info = gene + "_predicted_UTR\n"
    feature = "UTR"
    source = "TranStart"
    if direction == "+":
        UTR_start = str(TranStart)
        UTR_stop = str(start)
    if direction == "-":
        UTR_start = str(stop)
        UTR_stop = str(TranStop)
    new_gff_line = "\t".join([scaff, source, feature, UTR_start,
                             UTR_stop, score,
                             direction, frame, gene_info])
    return new_gff_line, UTR_start, UTR_stop


def avg_std_dev(positions):
    """function to return the avaerage and stadard deviation
    for a list of number.
    Uses Numpy to do the calculations"""
    # print("len pos = ", len(positions))
    # print(positions)
    if sum(positions) == 0:
        the_mean = 0
        standard_dev = 0
        return the_mean, standard_dev  
    try:
        the_mean = sum(positions) / float(len(positions))
        standard_dev = numpy.std(positions)
    except ValueError:
        the_mean = 0
        standard_dev = 0
    return the_mean, standard_dev


def mean_coverage(coverage_array, slice_start, slice_end):
    """function to get the mean coverage for a sliced region"""
    selected_coverage = coverage_array[slice_start : slice_end]
    return mean(selected_coverage)


def mean(list_of_values):
    """Calculate the mean average of a list of numbers."""
    # so don't have to worry about getting the divisor.
    # Explicit float(...) to allow for Python 2 division.
    try:
        mean = sum(list_of_values) / float(len(list_of_values))
        return mean
    except:
        return False

assert mean([1,2,3,4,5]) == 3


def middle_portion_of_transcript(seq):
    """return coordinates for middle portion"""
    coord_lower = int(len(seq)/4.0)
    coord_upper = int(len(seq)-coord_lower)
    return coord_lower, coord_upper


def subproces_func(cmd):
    """func to run subporcess"""
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
    return pipe


def fill_in_zero_cov(all_coverage, depth_file):
    """ function to fil in the zero values returned by samtools depth,
    or not returned by"""
    f_in = open(depth_file, "r")
    for line in f_in:
        # print(depth_filename)
        ref, possition, coverage = line.rstrip("\n").split("\t")
        possition = int(possition) - 1
        # assign the correct coverage to the relavent postiton,
        # thus removing the zero value.
        # Or if there is no info, it remains as a zero.
        all_coverage[possition] = int(coverage)
    f_in.close()
    return all_coverage


def run_samtools_depth(scaffold_start_stop, bam_file, outfile, logger):
    """funct to run samtools depth"""
    cmd = " ".join(["samtools",
                    "depth",
                    "-r",
                    scaffold_start_stop.rstrip(),
                    bam_file,
                    ">",
                    outfile])
    logger.info("samtools command: %s", cmd)
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
    return current_start, current_end


def add_one_direct_aware(current_start, current_end, interation_value,
                               direction):
    """function to alter the current start and end values based on the direction"""
    if direction == "+":
        current_start = current_start - interation_value
        current_end = current_end - interation_value
    else:
        current_start = current_start + interation_value
        current_end = current_end + interation_value
    return current_start, current_end


def iterate_coordinate_dict(gene_gff_line,
                            gene_name,
                            scaffold,
                            current_start,
                            current_stop,
                            logger):
    """check to see if the predicted TSS fall within a 
    predicted gene. It wanrs the user if the desired upstream
    regions falls into an existing gene. The coordinates will
    be altered upto this gene that is hits."""
    current_start = int(current_start)
    current_stop = int(current_stop)
    for gene, vals in gene_gff_line.items():
        # find the genes on the same scaffold
        if scaffold in vals:
            scaff, source, feature, start, stop, score, \
              direction, frame, gene_info = vals.split("\t")
            # if its is the same gene as the stop
            if gene_name in gene_info:
                # print(gene_name, gene_info)
                # we dont want to test it against intself!!
                continue
            if scaffold.rstrip() != scaff.rstrip():
                continue
            start = int(start) + 1
            stop = int(stop)
            # basically does the coordinate fall in the current
            # coordinate for a gene
            for UTR_coordinate in range(current_start, current_stop):
                ##logger.info(info)
                if UTR_coordinate > start and UTR_coordinate < stop:
                    warn = " ".join([gene_name,
                                     "UTR region falls into",
                                     "genic regions of",
                                     gene,
                                     "on scaffold",
                                     scaffold])
                    logger.warning(warn)
                    return "HITS genic region"
    return "OK"


## main function
def TranscriptionFind(genome, gene_start_stop_dict,
                      gene_first_exon_dict, gene_scaff_dict,
                      gene_direction, gene_set, gene_gff_line,
                      bam, stand_dev_threshold, walk, min_value,
                      interation_value, out_file, logger, TITLE,
                      keep_gene_depth,
                      default_cutoff,
                      test_mode):
    """iterate through gene in gff. Call samtools, get expression
    Does the expression fall within X standard deviations when compare
    to the first exon?
    """
    logger.info("RESULTS in outfile: %s", out_file)
    genome_index = index_genome_file(genome, logger)
    depth_set = set([])
    # open outfile:
    out_file_gff = out_file.split(".")[0] + "based_on_min_value.gff"
    out_file_gff2 = out_file.split(".")[0] + "based_on_SD_threshold.gff"
    logger.info("gff info will be in %s", out_file_gff)
    gff_out = open(out_file_gff, "w")
    gff_sd_out = open(out_file_gff2, "w")
    gene_failed_count = 0
    gene_results_printed_count = 0
    fall_off_contig_count = 0
    default_cutoff = float(default_cutoff)
    logger.info("any problem and default cutoff is used. Which is %.1f",
                default_cutoff)

    with open(out_file, 'w') as file_out:
        file_out.write(TITLE)
        for gene in gene_set:
            gene = gene.rstrip()
            start_stop = gene_start_stop_dict[gene]
            start, stop = start_stop.split("\t")
            start =int(start)
            stop = int(stop)
            scaffold = gene_scaff_dict[gene]
            scaffold = scaffold.rstrip()
            direction = gene_direction[gene]
            exon_start_exon_stop = gene_first_exon_dict[gene]
            exon_start, exon_stop = exon_start_exon_stop.split("\t")
            exon_start =int(exon_start)
            exon_stop = int(exon_stop)
            # call samtools to get the depth per posititon for
            # the transcript of interest
            depth_filename = os.path.join("temp_reads_per_base",
                                          gene + "_depth.tmp")
            #exon_depth_file = os.path.join("temp_reads_per_base",
                                          #gene + "_exon_depth.tmp")
            scaffold_depth_file = os.path.join("temp_reads_per_base",
                                               scaffold + "_depth.tmp")
            scaffold_start_stop = "%s:%s-%s" %(scaffold, start, stop)
            # call the func to run
            if scaffold_depth_file not in depth_set:
                depth_set.add(scaffold_depth_file)
                # print("not seen %s" % scaffold)
                pipe = run_samtools_depth(scaffold, bam_file,
                                          scaffold_depth_file, logger)
            # call the depth for the gene specifically
            pipe = run_samtools_depth(scaffold_start_stop, bam_file,
                                      depth_filename, logger)
            if "Y" not in keep_gene_depth.upper():
                # can keep the gene depth file, or not
                os.remove(depth_filename)

            # assign zeros to all positions of the transcript,
            # as samtool does no report zeros
            seq_record = genome_index[scaffold]
            if "Y" in test_mode.upper():
                logger.info("scaff = ", scaffold, "len scaffold = ", len(seq_record.seq),
                            "gene = ", gene, "depth scaff out = ", scaffold_depth_file)

            all_coverage = [0] * len(seq_record.seq)
            if "Y" in test_mode.upper():
                logger.info(" len all cov = ", len(all_coverage),
                            "all cov first 10 = ", all_coverage[:10])
            all_coverage = fill_in_zero_cov(all_coverage,
                                            scaffold_depth_file)
            # print("seq = ", len(seq_record.seq))
            # print(exon_all_coverage)
            # get the mean and std reads per base for exon 1
            exon_mean, exon_stdDev = avg_std_dev(all_coverage
                                                [exon_start:exon_stop])
            # get the mean and std reads per base for exon 1
            gene_mean, gene_stdDev = avg_std_dev(all_coverage
                                                [start:stop])
            if exon_mean == 0:
                warn = "No RNAseq expression for gene exon 1 %s" % gene
                logger.warning("%s: gene failed", warn)
                gene_failed_count = gene_failed_count + 1
                continue
            out_str = "\t".join([gene + ":",
                                "Cov min: %i" % min(all_coverage),
                                "max: %i" % max(all_coverage),
                                "gene mean %0.2f:" % gene_mean,
                                "gene std %0.2f:" % gene_stdDev,
                                "Sliced section:",
                                "exon mean %0.2f" % exon_mean,
                                "exon std: %0.2f" % exon_stdDev,
                                direction])
            # logger.info(out_str)
            cut_off = exon_mean - (int(stand_dev_threshold) * exon_stdDev)
            position_mean_cov = mean(all_coverage[exon_start:exon_stop])
            # walk in 3 bases to find the position where coverage sig drops
            current_end = stop
            current_start = start
            position_mean_cov = 10000000000000
            if cut_off < default_cutoff:
                logger.warning("%s gene cut off set to %.1f", gene, default_cutoff)
                cut_off = default_cutoff
            write = "yes"
            while position_mean_cov >= cut_off:
                current_start, current_end = walk_away_from_start(current_start,
                                                                  current_end,
                                                                  direction, walk)
                current_start, current_end = add_one_direct_aware(current_start,
                                                                  current_end,
                                                                  interation_value,
                                                                  direction)
                if current_start < 1:
                    logger.warning("%s has fallen off start scaffold %s",
                                   gene,
                                   scaffold)
                    position_mean_cov = 0
                    write = "no"
                    fall_off_contig_count = fall_off_contig_count + 1
                    break
                if current_end >= len(seq_record.seq):
                    logger.warning("%s has fallen off end scaffold %s",
                                   gene,
                                   scaffold)
                    position_mean_cov = 0
                    write = "no"
                    fall_off_contig_count = fall_off_contig_count + 1
                    break
                position_mean_cov = mean(all_coverage
                                         [current_start:current_end])
                if position_mean_cov == False:
                    position_mean_cov = 0
            #print("setting position_mean_cov to: ", position_mean_cov)
            # run the while loop again to find the position where the expression
            # is less than the option min value
            current_end1 = stop
            current_start1 = start
            position_mean_cov = 10000000000
            write_min_value = "ok"
            while position_mean_cov >= int(min_value):
                current_start1, current_end1 = walk_away_from_start(current_start1,
                                                                    current_end1,
                                                                    direction, walk)
                current_start1, current_end1 = add_one_direct_aware(current_start1,
                                                                    current_end1,
                                                                    interation_value,
                                                                    direction)
                if current_start < 1:
                    logger.warning("%s has fallen off start scaffold %s", gene, scaffold)
                    position_mean_cov = 0
                    write_min_value = "not_ok"
                    break
                if current_end >= len(seq_record.seq):
                    logger.warning("%s has fallen off end scaffold %s", gene, scaffold)
                    position_mean_cov = 0
                    write_min_value = "not_ok"
                    break
                # print("bases = ", all_coverage[current_start1:current_end1], "\n")
                position_mean_cov = mean(all_coverage[current_start1:current_end1])
                if position_mean_cov == False:
                    position_mean_cov = 0
                    break

            out_str = "\t".join([gene,
                                 str(current_start),
                                 str(current_end),
                                 str(seq_record.seq[current_start:current_end]),
                                 str(current_start1),
                                 str(current_end1),
                                 str(seq_record.seq[current_start1:current_end1]),
                                 "%s" % start,
                                 "%s" % stop,
                                 "%0.2f" % gene_mean,
                                 "%0.2f" % gene_stdDev,
                                 "%0.2f" % exon_mean,
                                 "%0.2f" % exon_stdDev,
                                 direction,
                                 "\n"])
            if current_start1 > 0 and current_end1 > 0 and current_start > 0 and current_end  > 0:
                if write == "yes" and write_min_value == "ok":
                    # print("writing: ", out_str)
                    file_out.write(out_str)
                    GENE_gff = gene_gff_line[gene]
                    # for the min value approach

                    new_gff_line1, UTR_start, UTR_stop = create_gff_line(GENE_gff, gene,
                                                   current_start1,
                                                   current_end1)
                    Min_val_Hits_geneic_or_not = iterate_coordinate_dict(gene_gff_line,
                                                             gene,
                                                             scaffold,
                                                             UTR_start,
                                                             UTR_stop,
                                                             logger)
                    if Min_val_Hits_geneic_or_not == "HITS genic region":
                        gene_failed_count = gene_failed_count + 1
                        continue
                    if Min_val_Hits_geneic_or_not == "OK":
                        gff_out.write(new_gff_line1)
                    # for the standard dev approach
                    new2_gff_line, UTR_start, UTR_stop = create_gff_line(GENE_gff, gene,
                                                    current_start,
                                                    current_end)
                    # Check to see if this hits a gene or not
                    sd_geneic_or_not = iterate_coordinate_dict(gene_gff_line,
                                                               gene,
                                                               scaffold,
                                                               UTR_start,
                                                               UTR_stop,
                                                               logger)
                    if sd_geneic_or_not == "HITS genic region":
                        gene_failed_count = gene_failed_count + 1
                        continue
                    if sd_geneic_or_not == "OK":
                        gff_sd_out.write(new2_gff_line)
                        gene_results_printed_count = gene_results_printed_count + 1
            else:
                gene_failed_count = gene_failed_count + 1
                
        logger.info("deleting scaffold depth files")
        for depthfile in depth_set:
            os.remove(depthfile)            
        logger.info("main function finished. %d gene failed", gene_failed_count)
        logger.info("Results generated for . %d gene", gene_results_printed_count)
        logger.info("fall_off_contig_count = %d", fall_off_contig_count)
        gff_out.close()


###############################################################################################
#to run it:


usage = """Use as follows:

Transcription start finder based on read depth coverage v0.1.0
THIS IS NOT finding the start
# of coding squence. but where transcription binds and starts

python TransStart.py -g genome.fasta
    --bam index_sorted_bam_file.bam
    --walk 5 --interation_value 1 
    --gff genes.gff -o outfile


Requirements:
    you must have in your PATH
    samtools
    Biopython
    numpy


steps:

1) Map RNAseq using STAR or similar tool
2) sort your bam out!
    samtools sort unsorted.bam sorted.bam
3) Index your sorted.bam
    samtools index sorted.bam


TranStart:
1 CPU time: 27, 000 genes ~ 5 hours
and requires 0.5 GB of RAM using a 500Mbp genome.

"""



parser = OptionParser(usage=usage)


parser.add_option("--gff", dest="gff",
                  default=None,
                  help="the predicted coordinate for the cds predictions .gff",
                  metavar="FILE")

parser.add_option("-g", "--genome",
                  dest="genome",
                  default=None,
                  help="the genome sequence. Not currently used. TO DO",
                  metavar="FILE")

parser.add_option("-b",
                  "--bam", dest="bam",
                  default=None,
                  help="the sorted, indexed bam file as a result of " +
                  "the reads being mapped back to the transcriptome  .bam",
                  metavar="FILE")

parser.add_option("-s",
                  "--std",
                  dest="stand_dev_threshold",
                  default=4,
                  help="If the expression of the current region is less then "
                  " this number of std away from the first exon mean. " +
                  "Default = 4")

parser.add_option("-w",
                  "--walk",
                  dest="walk",
                  default=3,
                  help="the number of bases to walk away to find expression " +
                  "Default = 3")

parser.add_option("--min_value",
                  dest="min_value",
                  default=2,
                  help="the min expression to return the coordinates for " +
                  "Default = 2")

parser.add_option("-i",
                  "--interation_value",
                  dest="interation_value",
                  default=1,
                  help="size of window to walk " +
                  "Default = 1")

parser.add_option("--default_cutoff",
                  dest="default_cutoff",
                  default=1,
                  help="if there is a problem with " +
                  "cut_off = exon_mean - (int(stand_dev_threshold) * exon_stdDev) ." +
                  "If this is less than default_cutoff, then " +
                  "default_cutoff value is applied. " +
                  " Default = 1")

parser.add_option("--help_full", dest="help_full",
                  default=None,
                  help="prints out a full description of this program")

parser.add_option("--keep_gene_depth",
                  dest="keep_gene_depth",
                  default="no",
                  help="keep the output of the depth for the genes" +
                  " yes or no ")

parser.add_option("--logger", dest="logger",
                  default=None,
                  help="Output logger filename. Default: " +
                  "outfile_std.log",
                  metavar="FILE")

parser.add_option("-o", "--out",
                  dest="outfile",
                  default="results.out",
                  help="Output filename (default: results.out)",
                  metavar="FILE")

parser.add_option("--test",
                  dest="test_mode",
                  default="No",
                  help="testing mode yes or no")

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
# --walk
walk = int(options.walk)
# --interation_value
interation_value = int(options.interation_value)
#--min_value
min_value = int(options.min_value)
# --default_cutoff
default_cutoff = options.default_cutoff
# test_mode
test_mode = options.test_mode

if not logger:
    log_out = "%s_%s_std.log" % (outfile, stand_dev_threshold)

if not os.path.isfile(bam_file):
    sys.exit("Input BAM file not found: %s" % bam_file)


#######################################################################
# Run as script
# Run as script
if __name__ == '__main__':
    # Set up logging
    if not options.logger:
        options.logger = "%s_%s_std.log" % (outfile, stand_dev_threshold)
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
    if not os.path.isfile(genome):
        logger.warning("Input genome file not found: %s" % genome)
        sys.exit("Input genome file not found: %s" % genome)
    if not os.path.isfile(bam_file):
        logger.warning("Input BAM file not found: %s" % bam_file)
        sys.exit("Input BAM file not found: %s" % bam_file)
    if os.path.isfile(bam_file + ".bai"):
        if os.path.getctime(bam_file + ".bai") < os.path.getctime(bam_file):
            logger.warning("your bam index file is older than bam")
            logger.warning("these must NOT match, so indexing again")
            cmd = " ".join(["samtools",
                            "index",
                            bam_file])
    if not os.path.isfile(bam_file + ".bai"):
        logger.warning("you have not indexed you bam file. I will do it.")
        cmd = " ".join(["samtools",
                        "index",
                        bam_file])
        # call the func
        pipe = subproces_func(cmd)
    # make a temp_folder_for_all_the_out_files
    dest_dir = os.path.join(os.getcwd(),
                            'temp_reads_per_base')
    try:
        os.makedirs(dest_dir)
    except OSError:
        print ("folder already exists")
    logger.info("Indexidng gff: %s", gff)
    gene_start_stop_dict, gene_first_exon_dict,\
                          gene_scaff_dict, gene_direction, gene_set,\
                          gene_gff_line = index_gff(gff, logger)
    logger.info("running analysis, running with the devil")
    TITLE = "\t".join(["#gene",
                       "TranStart(stats) Start",
                       "stop",
                       "sequence",
                       "TranStart(MinVal) Start",
                       "stop",
                       "sequence",
                       "gene_start",
                       "gene_stop",
                       "gene_mean_cov",
                       "gene_std_cov",
                       "exon mean",
                       "exon std",
                       "coding direction\n"])
    # main function
    TranscriptionFind(genome,
                      gene_start_stop_dict,
                      gene_first_exon_dict,
                      gene_scaff_dict,
                      gene_direction,
                      gene_set,
                      gene_gff_line,
                      bam_file,
                      stand_dev_threshold,
                      walk,
                      min_value,
                      interation_value,
                      outfile,
                      logger,
                      TITLE,
                      options.keep_gene_depth,
                      default_cutoff,
                      test_mode)
    out_file_gff = outfile.split(".")[0] + "based_on_min_value.gff"
    out_file_gff2 = outfile.split(".")[0] + "based_on_SD_threshold.gff"
    sort_cmd = "sort -k1n -k9n %s > temp1.gff" % out_file_gff
    sort_cmd2 = "sort -k1n -k9n %s > temp2.gff" % out_file_gff2
    mv_cmd = "mv temp1.gff %s" % out_file_gff
    mv_cmd2 = "mv temp2.gff %s" % out_file_gff2
    logger.info("sortintg gff output.")
    command_list = [sort_cmd, sort_cmd2, mv_cmd, mv_cmd2]
    for command in command_list:
        subproces_func(command)
    os.remove("example.log")
    logger.info("Pipeline complete: %s", time.asctime())

