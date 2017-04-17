#!/usr/bin/env python

############################################################################
# Title: possibly fix 5 prime ends of cds based of RNAseq coverage
############################################################################
"""
why? Sometime the 5 prime end is not correctly predcited.
This is really important to our
research. Can we imporve this?
"""
#############################################################################
# imports
import sys
import os
from optparse import OptionParser
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy
import subprocess
import tempfile
from collections import deque
import datetime
import logging
import logging.handlers
import time
# TODO make re look for Kozak
import re


if "-v" in sys.argv or "--version" in sys.argv:
    print("Fix five prime cds based on read depth coverage v0.1.0")
    cmd = "samtools 2>&1 | grep -i ^Version"
    sys.exit(os.system(cmd))

if "--help_full" in sys.argv:
    print("""Fix five prime cds based on read depth coverage v0.1.0:
why? Sometime the 5 prime end is not correctly predcited. This is really
important to our research. Can we imporve this?

Requires:
samtools 1.2 or greater
Biopython
numpy
some of these functions were taken from Peter Cock:

https://github.com/peterjc/pico_galaxy/tree/master/tools/coverage_stats

steps:

1) Index your transcriptome and map your reads back to your transcriptome:
    bowtie-build -f test_sequences.fasta test_sequences
    bowtie -S -p 8 test_sequences -1 R1.fq -2 R2.fq test_mapped.sam
2) Samtool index your transcriptome
    samtools faidx transcriptome.fasta
3) convert sam to bam
    samtools view -S -b -o unsorted.bam test_mapped.sam
4) sort your bam out!
    samtools sort unsorted.bam sorted.bam
5) Index your sorted.bam
    samtools index sorted.bam

How?

This script looks at the number of reads that map per base for your transcript
of interest.
If the number of mapped reads is greater than "threshold (default=100, reads
mapped per current base)" then it performs statistical
analysis on this to determine if the starting methionoine (as found in the
current predicted cds)
has significantly lower expression that the rest of the transcript. The
script will take the avergae coverage per base of this A,T,G
starting codon and compare this coverage to the mean of the middle 50%
of the transcript. If this sum(ATG)/3 has less then "threshold" number of
standard deviations of coverage less than the middle 50% of the transcript
then this may not be the correct starting codon.

As said, this "ATG" may not be the coreect starting postion. The program
will then looks for the next "ATG"
and performs statistical analysis to determine if this is a sensible starting
postiton. If the next starting ATG is within the threshold number of standard
deviations
from the mean of the middle 50% of the transcript, then this is set to the
new starting codon.

If this criteria is not met. Nothing is changed!

TO DO:
perform the same logic on the stop codon. Perform stats on the coverage
across the transcrip to identify fusions based on significantly different
expression profiles.
Travis testing. Codcov, pep8.

FIX five prime start in genomic regions based on this methodology

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
##############################################################################

try:
    # New in Python 3.4
    from statistics import mean
except ImportError:
    def mean(list_of_values):
        """Calculate the mean average of a list of numbers."""
        # so don't have to worry about getting the divisor.
        # Explicit float(...) to allow for Python 2 division.
        return sum(list_of_values) / float(len(list_of_values))

assert mean([1,2,3,4,5]) == 3


def sys_exit(msg, error_level=1):
   """Print error message to stdout and quit with given error level."""
   sys.stderr.write("%s\n" % msg)
   sys.exit(error_level)


def get_stats_no_window(depth_iterator, length):
    """Calculate min/max/mean coverage.

    Returns tuple (int, int, float).
    """
    # First call to iterator may return NoCoverage
    # This also makes initialising the min value easy
    try:
        ref, pos, depth = next(depth_iterator)
    except NoCoverage:
        return 0, 0, 0.0
    except StopIteration:
        raise ValueError("Internal error - was there no coverage?")
    total_cov = min_cov = max_cov = depth

    # Could check pos is strictly increasing and within 1 to length?
    for ref, pos, depth in depth_iterator:
        total_cov += depth
        min_cov = min(min_cov, depth)
        max_cov = max(max_cov, depth)
    mean_cov = total_cov / float(length)
    assert min_cov <= mean_cov <= max_cov
    return min_cov, max_cov, mean_cov


def get_stats_window(depth_iterator, length, window_size):
    """Calculate min/max/mean and min/max windowed mean.

    Assumes the depth_iterator will fill in all the implicit zero
    entries which ``samtools depth`` may omit!

    Assumes window_size < number of values in iterator!
    """
    window = deque()
    total_cov = 0
    min_cov = None
    max_cov = 0.0

    assert 1 <= window_size <= length

    prev_pos = 0
    while len(window) < window_size:
        try:
            ref, pos, depth = next(depth_iterator)
        except NoCoverage:
            return 0, 0, 0.0, 0.0, 0.0
        except StopIteration:
            outstr = "Not enough depth values to fill %i window" % window_size
            logger.info(outstr)
            raise ValueError("%s" % outstr)
        prev_pos += 1
        assert pos == prev_pos, "Discontinuity in cov vals for %s position %i" % (ref,
                                                                                  pos)
        total_cov += depth
        if min_cov is None:
            min_cov = depth
        else:
            min_cov = min(min_cov, depth)
        max_cov = max(max_cov, depth)
        window.append(depth)

    assert len(window) == window_size
    min_win = max_win = mean(window)
    for ref, pos, depth in depth_iterator:
        prev_pos += 1
        assert pos == prev_pos, "Discontinuity in cov val for %s position %i" % (ref,
                                                                                 pos)
        total_cov += depth
        min_cov = min(min_cov, depth)
        max_cov = max(max_cov, depth)
        window.popleft()
        window.append(depth)
        assert len(window) == window_size
        win_depth = mean(window)
        min_win = min(min_win, win_depth)
        max_win = max(max_win, win_depth)

    mean_cov = total_cov / float(length)

    assert prev_pos == length, "Missing final coverage?"
    assert len(window) == window_size
    assert min_cov <= mean_cov <= max_cov
    assert min_cov <= min_win <= max_win <= max_cov

    return min_cov, max_cov, mean_cov, min_win, max_win


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
    return_code = os.system(cmd)
    if return_code:
        clean_up()
        sys_exit("Return code %i from command:\n%s" % (return_code,
                                                       cmd))
    # creat a dictioanry to hold all the total expression values for the transcripts.
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


def strip_to_match_transcript_name(identifier):
    """remove the cds and split at the first pipe.
    This is to make the transcriptome and cds file names/ seq_record.id
    match"""
    if "cds." in identifier:
        # legacy format
        # e.g. >cds.Mp_O_10001_c0 etc ...
        return identifier.replace("cds.", "").split("|m.", 1)[0]
    else:
        # new transdecoder format.
        # e.g. >Mp_O_10001_c0::Mp_O_10001_c0_seq1::
        # we want the 2nd computer's [1] element
        return identifier.split("::")[1]


def find_longest_components(filename1,
                            cds_database,
                            out_filename):
    """this is a function to open up a fasta file, and
    producea a list of the longest representative transcripts per gene.
    This is only called is there are duplicated found with the same
    prefix name.
    Returns a new cds database with the longest cds per transcript only."""
    # this is out list of so called longest matches which we will
    # append and remove as applicable
    top_hits = []
    # current sequence length score value "to beat"
    current_length = int(0)
    # set up variables to assgn lastest values to ...
    transcriptome_Genes_names = set([])
    last_gene_name = ""
    last_component = ""
    loop_count = 0
    for seq_record in SeqIO.parse(filename1,
                                  "fasta"):
        sequence_len = len(seq_record)
        sequence_name = seq_record.id
        component = strip_to_match_transcript_name(sequence_name)
        # first time we see any record, save the values:
        if loop_count == 0:
            loop_count = loop_count + 1
            last_gene_name = sequence_name
            current_length = sequence_len
            last_component= component
            top_hits.append(seq_record.id)
        #########################################
        # first block: if the names are the same,
        # is the new length of sequence longer?
        if component == last_component:
            # print ("yes:", component, "component",  last_component,
            # "last_component", seq_record.id)
            # print ("current_length", current_length)
            if sequence_len > current_length:
                # print ("sequence_len > current_length", sequence_len,
                         #current_length)
                del top_hits[-1]
                top_hits.append(seq_record.id)
        ##########################################################
        # second block: if the name is new, put it in the name set.
        # use this sequence-length as the new one to "beat"
        else:
            top_hits.append(seq_record.id)
            last_gene_name = sequence_name
            current_length = sequence_len
            last_component= component
    outfile = open(out_filename, "w")
    for i in top_hits:
        seq_record = cds_database[i]
        SeqIO.write(seq_record, outfile,
                    "fasta")
    outfile.close()
    cds_database_new = SeqIO.index(out_filename,
                                   "fasta",
                                   key_function=strip_to_match_transcript_name)
    return cds_database_new


def parse_predicted_CDS_file(cds_file):
    """parse the cds file and index it.
    Take in the cds file and uses biopython to index it"""
    # this is for transdecoder names. May need to alter for other tools
    try:
        cds_database = SeqIO.index(cds_file,
                                   "fasta",
                                   key_function=strip_to_match_transcript_name)
        return cds_database
    except ValueError:
        outstr = ("WARNING: multi cds were predicted per transcript \n" +
        "\t- cannot change names. Going to pick the longest representative \n" +
        "\tcds per transcripts. Only do this if there are multiple cds \n" +
        "\tpredicted per transcript, otherwise message is not shown\n")
        logger.info(outstr)
        cds_database = SeqIO.index(cds_file, "fasta")
        # basically there are duplicates for each transcript.
        # So, find the longest representative and
        #  create a new cds_database, based on that
    # call function
    longest_rep = os.path.join("temp_fix_five_prime",
                               "longest_representative_seq.fasta")
    cds_database_new = find_longest_components(cds_file,
                                               cds_database,
                                               longest_rep)
    # return a seq_record object that can be
    # accessed in a dictionary like manner
    return cds_database_new


def index_genome_file(genome):
    """index the genome file
    return a seq_record object that can be accessed
    in a dictionary like manner"""
    genome_database = SeqIO.index(genome, "fasta")
    return genome_database


def index_gff_file(gff_file):
    """transcriptomes/cds prediction sometimes output GFF files
    this function indexes it. later these coordinates could be helpful.
    # return a dictionary. Key[transcript_name], vals are a list
    # containing
    # coordinates:  ['transdecoder', 'CDS', '1', '201', '.', '-', '.',
    # 'ID=cds.Mp_O_0_c0_seq1|m.2;Parent=Mp_O_0_c0_seq1|m.2'].
    Now looks like this:
    Mp_O_10020_c1_seq2 transdecoder gene 1 568	. + . I
                        D=Mp_O_10020_c1::Mp_O_10020_c1_seq2:
    """
    indexed_gff_file = dict()
    with open(gff_file, "r") as handle:
        # print ("i am here")
        for line in handle:
            if not line.strip():
                continue #if the last line is blank
            if line.startswith("#"):
                continue
            transcript_data = line.rstrip("\n").split()
            name = transcript_data[0]
            if transcript_data[2] == "CDS":
                # print "name = ",name
                indexed_gff_file[name]= transcript_data[1:]
    # return a dictionary. Key[transcript_name], vals are a list containing
    # coordinates:  ['transdecoder', 'CDS', '1', '201', '.', '-', '.',
    # 'ID=cds.Mp_O_0_c0_seq1|m.2;Parent=Mp_O_0_c0_seq1|m.2']
    return indexed_gff_file


def middle_portion_of_transcript(seq):
    """return coordinates for middle portion"""
    coord_lower = int(len(seq)/4.0)
    coord_upper = int(len(seq)-coord_lower)
    return coord_lower, coord_upper


def average_standard_dev(positions):
    """function to return the avaerage and stadard deviation
    for a list of number.
    Uses Numpy to do the calculations"""
    the_mean = sum(positions) / float(len(positions))
    standard_dev = numpy.std(positions)
    return the_mean, standard_dev


# function to walk along the transcript inframe as decided by the current
# reading frame for the cds
# and find the next "ATG (+)" or
def find_positive_next_ATG_b(transcriptome_record, position, strand):
    """function to find the next ATG.
    It is strand aware so will look in the correct direction.
    Returns the coordinate of the next methionine
    # function to walk along the transcript inframe as decided by the current
    # reading frame for the cds
    # and find the next "ATG (+)" or"""
    # print ("start position = %d" % position)
    transcriptome_record = transcriptome_record[position:]
    if strand == "+":
        next_codon = ""
        start = 0
        end = 3
        for i in range(len(transcriptome_record) - 60):
            start = start + 3
            end = end + 3
            outstr = ("start, end: ", start, end)
            logger.info(outstr)
            next_codon = transcriptome_record[start:end]
            # print next_codon
            # next_codon = transcript.seq[position+3:position+6]
            if next_codon == "ATG":
                next_methionine = start
                return next_methionine


def find_positive_next_ATG(transcriptome_record, position, strand):
    """function to find the next ATG; returns integer or None."""
    # print ("start position = %d" % position)
    transcriptome_record = transcriptome_record[position:]
    if strand =="+":
        next_codon = ""
        start = 0
        end = 3
        for i in range(len(transcriptome_record) - 60):
            start += 3
            end += 3
            next_codon = transcriptome_record[start:end]
            # print("start %i, end %i, codon: %s" % (start, end, next_codon))
            # print next_codon
            # next_codon = transcript.seq[position + 3:position + 6]
            if next_codon == "ATG":
                next_methionine = start
                return next_methionine + position  # Note restoring offet!
        return None  # To indicate no start codon found
        #raise ValueError("No start codon found")
    else:
        raise ValueError("Called on strand %r" % strand)


def find_downstream_start(name,
                          transcript,
                          current_start,
                          strand):
    """function to call other functions to find the next ATG
    start site.
    Takes in the transcript and the current
    start codon and strand coding direction"""
    if strand == "+":
        outstr = ("Looking for ATG after %d in seq: %s" % (current_start,
                                                           name))
        logger.info(outstr)
        if transcript[current_start:current_start+3] != "ATG":
            outstr = ("WARNING - existing annotation for " +
                      " %s does not start ATG" % name)
            logger.info(outstr)
        return find_positive_next_ATG(transcript, current_start, strand)
    elif strand == "-":
        outstr = ("Looking for CAT (i.e. ATG rev-comp) before " +
                  " %d in sequence %s" % (current_start, name))
        logger.info(outstr)
        new = find_positive_next_ATG(reverse_complement(transcript),
                                     len(transcript) - current_start, "+")
        if new is None:
            # No start codon found
            return None
        return len(transcript) - new
    else:
        raise ValueError("Bad strand value %r" % strand)


def mean_coverage(coverage_array, slice_start, slice_end):
    """function to get the mean coverage for a sliced region"""
    selected_coverage = coverage_array[slice_start : slice_end]
    return mean(selected_coverage)


def translate_cds(fasta_file):
    """function to translate the new cds file"""
    # outfile is define in optparser
    file_out_name = outfile + ".pep"
    file_out = open(file_out_name, "w")
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        seq_record.seq = seq_record.seq.translate()
        SeqIO.write(seq_record, file_out, "fasta")
    file_out.close()


##### main function
def parse_transcriptome_file(genome,
                             transcriptome_file,
                             cds_file,
                             bam,
                             gff,
                             min_read_count,
                             min_max_cov_per_base,
                             stand_dev_threshold,
                             outfile,
                             overall_expression_dic,
                             out_file):
    """iterate through transcriptome. Call samtools, get expression
    Does the expression fall within X standard deviations when compare
    to the middle 50% of the transcript?
    """
    # call the function to return a biopython seqIO idexed database
    # format the CDS file
    # note these names should have been altered to return names
    # in the same manner as the transcriptome
    cds_index_database = parse_predicted_CDS_file(cds_file)

    # open outfile:
    file_out = open(out_file, "w")

    # call function to get the overall transcript expression,
    # based on number of reads mapped:
    # calls samtools
    overall_expression_dic = get_total_coverage(bam,
                                                "overall_transcript" +
                                                "_expression.txt")

    # threshold for min_read_count
    min_read_count_threshold = int(min_read_count) # default is 30

    for transcriptome_record in SeqIO.parse(transcriptome_file,
                                            "fasta"):
        # get the overall expression from dict for the
        # transcriptome reocrd of interest
        transcript_expression = overall_expression_dic[transcriptome_record.id]
        # split up the line, basically this is how
        # samtools idxstats spits it out.
        # ID, Length, read count, last_coloumn
        length, read_count, last_coloumn = transcript_expression
        # santity check that the length match up
        assert length == len(transcriptome_record.seq)

        # basically if it has less than ~30? reads mapped
        # we cant really do stats on it
        if read_count < min_read_count_threshold:
            continue

        # test one: is the a cds predcited for this, if not move on
        try:
            cds_record = cds_index_database[transcriptome_record.id]
            original_cds_record = cds_index_database[transcriptome_record.id]
        except KeyError:
            # transdecoder uses cdhit 90% so some wont have cds...
            outstr = ("no CDS was found for " +
                      "%s . Moving to the next" % (transcriptome_record.id))
            logger.info(outstr)
            continue
        # call samtools to get the depth per posititon for
        # the transcript of interest
        depth_filename = os.path.join("temp_fix_five_prime",
                                      "depth.tmp")
        cmd = " ".join(["samtools",
                        "depth",
                        "-r",
                        transcriptome_record.id,
                        bam_file,
                        ">",
                        depth_filename])
        outstr = ("Running %s" % cmd)
        logger.info(outstr)
        return_code = os.system(cmd)
        assert return_code == 0, """samtools says NO!!
        - something went wrong. Is your BAM file correct?
        Samtools version? Did you sort you .bam?"""

        # assign zeros to all positions of the transcript,
        # as samtool does no report zeros
        all_coverage = [0] * len(transcriptome_record.seq)

        for line in open(depth_filename):
            ref, possition, coverage = line.rstrip("\n").split("\t")
            possition = int(possition) - 1

            # assign the correct coverage to the relavent postiton,
            # thus removing the zero value.
            # Or if there is no info, it remains as a zero.
            all_coverage[possition] = int(coverage)
        # Now have all the coverage information for this transcript in all_coverage
        if max(all_coverage)< int(min_max_cov_per_base):
            continue
        start_position = transcriptome_record.seq.find(cds_record.seq)

        if start_position == -1:
            # reverse complement the transcript for negative coding CDS
            transcriptome_record.seq = transcriptome_record.seq.reverse_complement()
            start_position = transcriptome_record.seq.find(cds_record.seq)

        end_position = start_position + len(cds_record.seq)
        assert 0 <= start_position < end_position <= len(transcriptome_record)
        coord_lower, coord_upper \
                = middle_portion_of_transcript(transcriptome_record.seq)
        outstr = ("coord_lower = " +
                  "%d\tstart_position = %d" %(coord_lower,
                                              start_position))
        logger.info(outstr)
        if coord_lower <= start_position:
            coord_lower = start_position+50
        if coord_lower + 30 >= coord_upper:
            continue

        the_mean, standard_dev = average_standard_dev (all_coverage
                                                       [coord_lower:coord_upper])
        out_str = " ".join([transcriptome_record.id + ":",
                            "Cov min: %i" % min(all_coverage),
                            "max: %i" % max(all_coverage),
                            "Sliced section: ",
                            "mean %0.2f" % the_mean,
                            "std: %0.2f" % standard_dev,
                            "(+) coding strand"])

        # to make it neater, I split the sentance over 2 strings
        logger.info(out_str)

        cut_off = the_mean - (int(stand_dev_threshold)
                              * standard_dev)
        #if all_coverage[start_position] < cut_off:
        start_codon_mean_cov = mean(all_coverage
                                   [start_position : start_position + 3])
        if start_codon_mean_cov < cut_off:
            out_str = " ".join(["houston we have a problem!!: ",
                                transcriptome_record.id,
                                " diff expression.",
                                "Cov min: %i" % min(all_coverage),
                                "max: %i" % max(all_coverage),
                                "For sliced section: mean: %0.2f" % the_mean,
                                "std: %0.2f" % standard_dev,
                                "current start position: %0.2f" % (all_coverage[start_position])])
            logger.info(out_str)

            #find the next ATG:
            next_ATG_position = find_downstream_start(transcriptome_record.id,
                                                      str(transcriptome_record.seq),
                                                      start_position, "+")
            # if it cant find another ATG - dont change the data....
            if next_ATG_position == None:
                # No alternative start, do nothing
                #reset record
                cds_record = original_cds_record
                pass
            elif next_ATG_position > end_position-60:
                # No alternative start within CDS, do nothing
                #reset record
                cds_record = original_cds_record
                pass
            else:
                assert next_ATG_position < end_position, "Was %i to %i, new start %i" % (start_position,
                                                                                         end_position,
                                                                                         next_ATG_position)
                if all_coverage[next_ATG_position] > cut_off:
                    out_str =  ("next_ATG_position  = %d" % (next_ATG_position))
                    logger.info(out_str)
                    out_str = " ".join(["houston, we may have fixed the problem!! ",
                                        transcriptome_record.id,
                                        "length: %d" % len(transcriptome_record.seq),
                                        "min %i, max %i, " % (min(all_coverage),
                                                              max(all_coverage)),
                                        "For sliced section: mean: %0.2f" % the_mean,
                                        "std: %0.2f" % standard_dev,
                                        "current start position: %0.2f" % (all_coverage[start_position]),
                                        "\n\t",
                                        "NEW start %i, end %i" % (next_ATG_position, end_position),
                                        "For sliced section:",
                                        "mean %0.2f, std %0.2f " % (the_mean, standard_dev),
                                        "next_ATG_position %0.2f" % (all_coverage[next_ATG_position])])
                    logger.info(out_str)
                    cds_record.seq = transcriptome_record.seq[next_ATG_position:end_position] # <-- TESTS work.
                    cds_record.description = "Five_prime_altered new coordinates: %d - %d\t" % \
                                             (next_ATG_position, end_position) \
                                             + cds_record.description
                    if len(cds_record.seq) < end_position - 60:
                            #reset to original
                            out_str = ("houston, im not resetting it!!!:")
                            logger.info(out_str)
                            out_str = " ".join(["the length of the",
                                                "cds_record would be",
                                                "%d," % len(cds_record.seq),
                                                "last acceptable coordinate ",
                                                "%d" % (end_position - 60)])
                            logger.info(out_str)
                            
                            cds_record = original_cds_record
                else:
                    another_ATG_position = find_downstream_start(transcriptome_record.id,
                                                                 str(transcriptome_record.seq),
                                                                 next_ATG_position, "+")
                    #if this new ATG comed witht the thresholds for "being" real
                    if another_ATG_position is None:
                        #reset to original
                        cds_record = original_cds_record
                        out_str = ("No alt start, cannot do anything with this %s" % (transcriptome_record.id))
                        logger.info(out_str)
                    elif another_ATG_position > end_position:
                        # No alternative start within CDS, do nothing
                        #reset to original
                        cds_record = original_cds_record
                        pass
                    elif all_coverage[another_ATG_position] > cut_off:
                        #set the new cds
                        cds_record.description = original_cds_record.description
                        cds_record.seq = transcriptome_record.seq[another_ATG_position:end_position]
                        out_str = " ".join(["houston, we may have fixed the problem!!2 ",
                                            transcriptome_record.id,
                                            "another_ATG_position_NEW start",
                                            "length: %d" % len(transcriptome_record.seq),
                                            "min %i, max %i, " % (min(all_coverage),
                                                                  max(all_coverage)),
                                            "For sliced section: mean: %0.2f" % the_mean,
                                            "std: %0.2f" % standard_dev,
                                            "current start position: %0.2f" % (all_coverage[start_position]),
                                            "\n\t",
                                            "NEW start %i, end %i" % (next_ATG_position, end_position),
                                            "For sliced section:",
                                            "mean %0.2f, std %0.2f " % (the_mean, standard_dev),
                                            "another_ATG_position %0.2f" % (all_coverage[next_ATG_position])])
                        logger.info(out_str)

                        cds_record.description = "Five_prime_altered2 new coordinates: %d - %d\t" % \
                                             (another_ATG_position, end_position) \
                                             + cds_record.description
                        if len(cds_record) > end_position - 60:
                            #reset to original
                            out_str = ("houston, im not resetting it!!!2:")
                            logger.info(out_str)
                            cds_record = original_cds_record

                    else:
                        out_str = ("cannot do anything %s -- leave cds as is" % (transcriptome_record.id))
                        logger.info(out_str)
                        #reset to original
                        cds_record = original_cds_record
        if not len(cds_record):
            #reset to original
            cds_record = original_cds_record

        assert len(cds_record), "Trimmed to nothing? %s " %(transcriptome_record.id)
        # out_str = ("im writing %s" %(transcriptome_record.id))
        # logger.info(out_str)
        SeqIO.write(cds_record, file_out, "fasta")


###############################################################################################
#to run it:


usage = """Use as follows:

python Fix_five_prime_CDS.py -t trnascriptome --cds nt_coding_seq --bam index_sorted_bam_file.bam
    --gff if_you_have_one --exp outfile_name_for_expression_values (default: overall_reads_mapped_per_sequences.txt)
    --prot (protein_cds_not_currently_used) -o outfile


Requirements:
    you must have samtools in your PATH
    Biopython
    numpy


why? Sometime the 5 prime end is not correctly predcited. This is really important to our
research. Can we imporve this?


some of these functions were taken from Peter Cock:

    https://github.com/peterjc/pico_galaxy/tree/master/tools/coverage_stats

steps:

1) Index your transcriptome and map your reads back to your transcriptome:
    bowtie-build -f test_sequences.fasta test_sequences
    bowtie -S -p 8 test_sequences -1 R1.fq -2 R2.fq test_mapped.sam
2) Samtool index your transcriptome
    samtools faidx transcriptome.fasta
3) convert sam to bam
    samtools view -S -b -o unsorted.bam test_mapped.sam
4) sort your bam out!
    samtools sort unsorted.bam sorted.bam
5) Index your sorted.bam
    samtools index sorted.bam


"""



parser = OptionParser(usage=usage)

parser.add_option("-t", "--transcriptome", dest="transcriptome", default=None,
                  help="the transcriptome assembly .fasta",
                  metavar="FILE")

#cds = options.cds
#default_gff = cds.replace("cds", "gff")
#default_pro = cds.replace("cds", "pep")

parser.add_option("--cds", dest="cds",
                  default=None,
                  help="the predicted cds from the transcriptome assembly .fasta",
                  metavar="FILE")

parser.add_option("--gff", dest="gff",
                  default=None,
                  help="the predicted coordinate for the cds predictions .gff",
                  metavar="FILE")

parser.add_option("--prot", dest="prot",
                  default=None,
                  help="the predicted amino acid cds  .pep",
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
parser.add_option("--min_read_count",
                  dest="min_read_count",
                  default="100",
                  help="the min_read_count that a transcript must have before it is"
                  "considered for statistical analysis. Default = 100")

parser.add_option("--min_max_cov_per_base",
                  dest="min_max_cov_per_base",
                  default="30",
                  help="the min_max_cov_per_base that a transcript must have before it is"
                  "considered for statistical analysis. Default = 30")

parser.add_option("--std",
                  dest="stand_dev_threshold",
                  default="3",
                  help="If the expression of the start of the cds is stand_dev_threshold"
                  " +/- the mean the flag as to be looked at. Default = 3")

parser.add_option("--help_full", dest="help_full", default=None,
                  help="prints out a full description of this program")

parser.add_option("--exp", dest="over_all_expression",
                  default="overall_reads_mapped_per_sequences.txt",
                  help="Output filename for the overall number of reads"
                  "that map to the sequence of interest. "
                  "Default: overall_reads_mapped_per_sequences.txt",
                  metavar="FILE")

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
#-t
transcriptome_file = options.transcriptome
#--cds
cds_file = options.cds
#--gff
gff = options.gff
#--bam
bam_file = options.bam
#-o
outfile= options.outfile
# -min_read_count
min_read_count = options.min_read_count
# ==stand_dev_threshold
stand_dev_threshold = options.stand_dev_threshold
# --exp
over_all_expression = options.over_all_expression
# --min_max_cov_per_base
min_max_cov_per_base = options.min_max_cov_per_base
# logfile
logger = options.logger
if not logger:
    log_out = "%s_%s_std.log" % (outfile, stand_dev_threshold)

# simple test
assert mean([1,2,3,4,5]) == 3

if not os.path.isfile(transcriptome_file):
    sys_exit("Input transcriptome file not found: %s" % transcriptome)

if not os.path.isfile(bam_file):
    sys_exit("Input BAM file not found: %s" % bam)


#######################################################################
# Run as script
if __name__ == '__main__':
    if not os.path.isfile(transcriptome_file):
        sys_exit("Input transcriptome file not found: %s" % transcriptome)
    if not os.path.isfile(cds_file):
        sys_exit("Input cds_file file not found: %s" % cds_file)
    if not os.path.isfile(bam_file):
        sys_exit("Input BAM file not found: %s" % bam)

    # Set up logging
    logger = logging.getLogger('Fix_five_prime_CDS.py: %s' % time.asctime())
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
        outstr = "Could not open %s for logging" % logger
        logger.error(outstr)
        sys.exit(1)
    # Report input arguments
    outstr = "Command-line: %s" % ' '.join(sys.argv)
    logger.info(outstr)

    overall_expression_dic = get_total_coverage(bam_file,
                                                over_all_expression)

    parse_transcriptome_file(genome,
                             transcriptome_file,
                             cds_file,
                             bam_file,
                             gff,
                             min_read_count,
                             min_max_cov_per_base,
                             stand_dev_threshold,
                             outfile,
                             overall_expression_dic,
                             outfile)
    print ("translating the new cds file")
    translate_cds(outfile)

