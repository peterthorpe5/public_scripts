#!/usr/bin/env python3
# title: Parse clusters and find the database species
# in the cluster
# author: Peter Thorpe and Leighton Pritchard
# September 2016. The James Hutton Insitute, Dundee, UK.

# imports
import os
import sys
from optparse import OptionParser
import datetime
import os
from sys import stdin,argv
import subprocess

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator


###############################################################

def stop_err( msg ):
    """stop function"""
    sys.stderr.write("%s\n" % msg)
    sys.exit(1)

try:
    # New in Python 3.4
    from statistics import mean
except ImportError:
    def mean(list_of_values):
        """Calculate the mean average of a list of numbers."""
        # Quick and dirty, assumes already a list not an
        # interator.
        # so don't have to worry about getting the divisor.
        # Explicit float(...) to allow for Python 2 division.
        return sum(list_of_values) / float(len(list_of_values))

assert mean([1,2,3,4,5]) == 3


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


def get_fasta_stats(fasta):
    """function to get stats on a given fasta file"""
    with open(fasta, 'r') as seq:
        sizes = [len(record) for record in \
                 SeqIO.parse(seq, 'fasta')]
    min_contig = min(sizes)
    max_contig = max(sizes)
    avg_contig = mean(sizes)
    num_contig = len(sizes)
    return min_contig, max_contig, avg_contig, num_contig


def count_reads(in_name, current_read_count):
    """function to count the number of reads found in a cluster.
    This is determined by the abundance count at the end of the read
    name """
    reads = in_name.split("_")[-1]
    try:
        reads = reads.split("abundance=")[1]
    except ValueError:
        # dont brake this
        # this is a database entry
        reads = 1
    return current_read_count + int(reads)


def parse_line(line):
    """finction to parse line"""
    if not line.strip():
        return False #if the last line is blank
    if line.startswith("#"):
        return False
    line = line.rstrip()
    out_put_str = ""
    if "\t" in line:
        cluster_line = line.rstrip("\n").split("\t")
    else:
        cluster_line = line.rstrip("\n").split()
    # are the database Phy singletons? - then we wont be interested
    number_of_reads_hitting_species = len(cluster_line)
    # set up some blank variables
    species = ""
    # database counter
    db_species = 0
    return cluster_line, number_of_reads_hitting_species,\
           species, db_species


def align_cluster(infile):
    """function to muscle align a cluster"""
    cmd = 'muscle -in %s -out %s_TEMP' % (infile, infile)
    print("Running %s" % cmd)
    return_code = os.system(cmd)
    assert return_code == 0, """muscle says NO!! -
                             is muscle in your PATH?"""
    # wait  forces Linux to wait until all command are complete
    cmd_wait = 'wait'
    print("Running %s" % cmd_wait)
    return_code = os.system(cmd_wait)
    # refine take longer but is the most accurate
    cmd2 = 'muscle -in %s_TEMP -out %s_aligned.fasta -refine' % (infile,
                                              infile.split(".fasta")[0])
    aligned_file_name = "%s_aligned.fasta" % (infile.split(".fasta")[0])
    print("Running %s" % cmd2)
    return_code = os.system(cmd2)
    assert return_code == 0, """muscle -refine says NO!! -
                                is muscle in your PATH, does
                                this version allow this parameter?"""
    cmd_wait = 'wait'
    print("Running %s" % cmd_wait)
    return_code = os.system(cmd_wait)

    cmd3 = 'rm %s*_TEMP' % (infile)
    print("Running %s" % cmd3)
    return_code = os.system(cmd3)
    assert return_code == 0, """something went wrond deleting
                                the temp alignment file!"""
    # passing the aligned file name back for seq_identity func
    return aligned_file_name


def get_names_from_Seq_db(seq_db):
    """function to get a list of name in the seq db"""
    names = []
    names_abudance_removed = []
    db = open(seq_db, "r")
    for seq_record in SeqIO.parse(db, "fasta"):
        if seq_record.id.endswith("_1"):
            names.append(seq_record.id)
            names_abudance_removed.append(("_").join(
                seq_record.id.split("_")[:-1]))
        else:
            names_abudance_removed.append(seq_record.id)
            names.append(seq_record.id + "_1")
    db.close()
    # print (names_abudance_removed)
    return names, names_abudance_removed


def coded_name_to_species(database_file):
    """functiong takes the already generated tab separated
    database of coded name to species file. Returns a dic
    of coded_name to species"""
    with open(database_file) as fh:
        data = fh.read().split("\n")
    coded_name_to_species_dict = dict()
    coded_species_to_name_dict = dict()
    for line in data:
        if not line.strip():
            continue  # if the last line is blank
        if line.startswith("#"):
            continue
        data = line.split("\t")
        coded_name = data[0]
        species = data[1:]
        coded_name_to_species_dict[coded_name.rstrip()] = species
        for i in species:
            coded_species_to_name_dict[i] = coded_name
    # print (coded_name_to_species_dict)
    return coded_name_to_species_dict, coded_species_to_name_dict


def write_out_cluster_as_fa(member, all_sequences,
                            all_sequences_no_abundance,
                            names,
                            names_abudance_removed,
                            coded_name_to_read_dict,
                            coded_species_to_name_dict,
                            cluster_fasta):
    """function to write out the cluster as
    a fasta file"""
    member = member.rstrip()
    db_name = ("_").join(member.split("_")[:-1])
    abundance = member.split("_")[-1]
    if member in names or member in names_abudance_removed or db_name in names_abudance_removed:
        # this is db entry
        try:
            seq_record = all_sequences[db_name]
        except:
            seq_record = all_sequences_no_abundance[db_name]
        SeqIO.write(seq_record, cluster_fasta, "fasta")
        return True
    else:
        hex_name = coded_species_to_name_dict[db_name]
        #print ("mem= ", member, "db=", db_name, "hex =", hex_name)
        #print (coded_species_to_name_dict[db_name])
        try:
            seq_record = all_sequences[hex_name]
        except:
            seq_record = all_sequences_no_abundance[hex_name]
        seq_record.description = ""
        seq_record.id = db_name + "_abunance=" + abundance
    if not seq_record:
        # break the program here
        print (coded_species_to_name_dict)
        stop_err ("""missing sequence %s.
                    This should not be seen.
                    Have you passed me the correct
                    fasta file?\n\n""" % member)
    try:
        SeqIO.write(seq_record, cluster_fasta, "fasta")
    except ValueError:
        print ("could not write record for %s" % member)
    return True


def strip_to_match_db_name(identifier):
    """function to remove the abundance off the names, so
    the keys in the dic will match later"""
    return ("_").join(identifier.split("_")[:-1])


# main function
def parse_tab_file_get_clusters(fasta, seq_db, all_fasta,
                                in_file, old_to_new,
                                min_novel_cluster_threshold,
                                right_total_reads,
                                working_directory, v,
                                align, out_file):
    """script to open up a tab separeted clustering output
    and identify the species in the clustering"""

    if seq_db:
        names, names_abudance_removed = get_names_from_Seq_db(seq_db)
    else:
        names = []

    # index the all_seq fasta. remove the abundance
    all_sequences_no_abundance =  SeqIO.index(all_fasta,
                                              "fasta",
                                              key_function=strip_to_match_db_name)
    all_sequences =  SeqIO.index(all_fasta, "fasta")

    # call function to get old to new name:
    coded_name_to_read_dict, \
        coded_species_to_name_dict = coded_name_to_species(old_to_new)
    # call the function to get fasta stats
    min_contig, max_contig, avg_contig, \
                num_contig = get_fasta_stats(fasta)
    cluster_file = open (in_file, "r")
    summary_out_file = open(out_file, "w")

    ITS_hitting_db_species = 0
    # title for the results file
    title = "# cluster_number\tspecies\t" +\
            "number_of_reads_hitting_species\n"
    summary_out_file.write(title)
    cluster_number = 0
    novel_count = 0
    for line in cluster_file:
        cluster_number += 1
        # call function to parse line
        if not parse_line(line):
            continue
        reads_of_interest = ""
        cluster_line, number_of_reads_hitting_species, \
                      species, db_species = parse_line(line)
        # print ("cluster_number = ", cluster_number,
               # "number_of_reads_hitting_species = ",
               # number_of_reads_hitting_species)
        # basically not a singleton cluster
        if len(cluster_line) > 1:
            # open a file to put seq into
            out_name = os.path.join(working_directory,
                                    "clusters_d%s" % str(v),
                                    "cluster%d_len%d.fasta" % (cluster_number,
                                                               len(cluster_line)))
            cluster_fasta = open(out_name, "w")
        for member in cluster_line:
            # write out the cluster members to a fasta file
            if len(cluster_line) > 1:
                write_out_cluster_as_fa(member, all_sequences,
                                        all_sequences_no_abundance,
                                        names,
                                        names_abudance_removed,
                                        coded_name_to_read_dict,
                                        coded_species_to_name_dict,
                                        cluster_fasta)

            # is a memeber of the database in this cluster?
            temp_name = ("_").join(member.split("_")[:-1])
            if temp_name in names or temp_name in names_abudance_removed:
                # yes we are interested in this cluster
                db_species = db_species + 1
                species = species + member + " "
                # remove all the memebr which are database memebers
                number_of_reads_hitting_species = number_of_reads_hitting_species - 1

        # another loop
        if number_of_reads_hitting_species >= 1 and db_species > 0:
            # after the line, are there more memeber that are
            # not database members?
            for member in cluster_line:
                if ("_").join(member.split("_")[:-1]) in names or member in names_abudance_removed:
                    continue
                # we do not add the reads that hit the database
                reads_of_interest = reads_of_interest + member + " "

            ITS_hitting_db_species = ITS_hitting_db_species + \
                                     number_of_reads_hitting_species
            # format the data

            data_output = "%d\t%s\t%d\n" %(cluster_number, species.rstrip(), \
                                           number_of_reads_hitting_species)

            # write out the data
            summary_out_file.write(data_output)
            db_reads_perc = (float(ITS_hitting_db_species)/ num_contig)*100
        else:
            novel_count = 0
            if len(cluster_line) > min_novel_cluster_threshold:
                # open a file to put seq into
                if db_species < 1:
                    novel_count = 1
                    out_novel = os.path.join(working_directory,
                                             "novel_d%s" % str(v),
                                             "novel%d_len%d.fasta" % (cluster_number,
                                                                      len(cluster_line)))

                    novel_fasta = open(out_novel, "w")
                    for member in cluster_line:
                        write_out_cluster_as_fa(member, all_sequences,
                                                all_sequences_no_abundance,
                                                names,
                                                names_abudance_removed,
                                                coded_name_to_read_dict,
                                                coded_species_to_name_dict,
                                                novel_fasta)
                else:
                    novel_fasta = False
            else:
                novel_count = 0

        try:
            cluster_fasta.close()
            if align:
                # align the cluster
                # print ("I am about to align the cluster: %s" % out_name)
                aligned_file_name = align_cluster(out_name)
        except ValueError: # singleton cluster - no file generated
            pass
        try:
            if novel_count > 0:
                novel_fasta.close()
            if align:
                # align the cluster
                # print ("I am about to align the cluster: %s" % out_novel)
                aligned_file_name = align_cluster(out_novel)
        except ValueError:  # no novel files to close
            pass
    db_reads_perc = (float(ITS_hitting_db_species)/ num_contig)*100

    fasta_file_summary = """    # Fasta file assembled seq summary:
    # min_contig = %d max_contig = %d avg_contig = %d
    # Total number of assembled sequences = %d
    # number of assembled-reads clustering with database = %d
    # number of starting reads = %d
    # percent of assembled-reads clustering with database = %.2f\n""" %(min_contig,
                                                    max_contig, avg_contig,
                                                    num_contig,
                                                    ITS_hitting_db_species,
                                                    right_total_reads,
                                                    db_reads_perc)
    summary_out_file.write(fasta_file_summary)

    dbmin_contig, dbmax_contig, dbavg_contig, \
                dbnum_contig = get_fasta_stats(seq_db)
    db_fa_summary = """    # db_fa_summary:
    # dbmin_contig = %d dbmax_contig = %d dbavg_contig = %d
    # number of db_seq = %d
    """ %(dbmin_contig, dbmax_contig, dbavg_contig,
          dbnum_contig)
    summary_out_file.write(db_fa_summary)

    cluster_file.close()
    summary_out_file.close()

#############################################################################
# to run the script

usage = """usage :

this scripts to summarise what clusters with what Phy species. The script also
dumps one file per cluster, BLASTN searches these back against itself
and aligns it.

WARNING: THIS WILL PRODUCE 10 000S OF THOUSANDS OF FILES WHICH WILL SLOW
YOU COMPUTER DOWN, AND MAKE YOU UNPOPULAR WITH SYS ADMIN. YOU ASKED FOR
THIS FEATURE!!! The BLAST and alignment also makes it very slow.


$ python parse_clusters.py -i clustering_outfile_already_decoded_from_temp_names
-n 3 -a all_seqeucnes.fasta -o summarise_clusters.out

command line option

requires Biopython

"""

parser = OptionParser(usage=usage)

parser.add_option("-f","--fasta", dest="fasta",
                  default=None,
                  help="assembled fasta file to get stats on",
                  metavar="FILE")

parser.add_option("-a","--all_fasta", dest="all_fasta",
                  default=None,
                  help="both the assembled fasta and databse fasta."
                  " this is to get the sequences from from for "
                  " the clusters folder",
                  metavar="FILE")

parser.add_option("--seq_db", dest="seq_db", default=False,
                  help="the databse used to cluster with")

parser.add_option("-n","--min_novel_cluster_threshold",
                  dest="min_novel_cluster_threshold", default=4,
                  help="min_novel_cluster_threshold to "
                  " determine is this is a "
                  " novel clusters to output. Default = 4 ")

parser.add_option("-l", "--left", dest="left",
                  default="temp_not_trimmedr1.fq",
                  help="left reads unzipped fq file. needed to "
                  " get count of reads. ",
                  metavar="FILE")

parser.add_option("--Name_of_project", dest="Name_of_project",
                  default="RESULTS.out",
                  help="name of project to make folders. ")

parser.add_option("-i","--in", dest="in_file", default=None,
                  help="clustering out file")

parser.add_option("-v","--difference", dest="v", default="1",
                  help="current swarm v value")

parser.add_option("-r", "--right", dest="right",
                  default="temp_not_trimmedr2.fq",
                  help="right reads unzipped fq file. needed to "
                  " get count of reads. ",
                  metavar="FILE")

parser.add_option("--old_to_new", dest="old_to_new", default=None,
                  help="file with the old and new names tab separated.",
                  metavar="FILE")

parser.add_option("--align", dest="align",
                  default=False,
                  help="this option performs alignment on a cluster. "
                  " Turn off for faster results."
                  " To turn off --align False"
                  " muscle must be in your PATH called muscle")

parser.add_option("-o", "--out_prefix", dest="out_file",
                  default="summarise_clusters.out",
                  help="prefix to the output filenames")


(options, args) = parser.parse_args()
# -f
fasta = options.fasta
# --all_fasta
all_fasta = options.all_fasta
# --seq_db
seq_db = options.seq_db
# -i
in_file = options.in_file
# -l
left = options.left
# -r
right = options.right
# -n
min_novel_cluster_threshold = options.min_novel_cluster_threshold
# -o
out_file = options.out_file
# --align
align = options.align
# -- old_to_new
old_to_new = options.old_to_new
# -v
v = options.v
# --Name_of_project
Name_of_project = options.Name_of_project


###################################################################
# start of program

# call the function to count the trimmed reads
right_total_reads = count_fq_reads(right)
left_total_reads = count_fq_reads(left)

# sanity test.
assert right_total_reads == left_total_reads, """ \n\nthe total
number of reads in the
left and right file do not match.
Something has gone wrong before here! These should
be the same. Are you giving me the the correct files?\n\n"""

#fasta, seq_db, all_fasta,
# make sure this is an int
min_novel_cluster_threshold = int(min_novel_cluster_threshold)
v = int(v)

###################################################################
# make folders to put the cluster fasta files.
make_folder_list = ["clusters", "novel"]
working_dir = os.getcwd()
for name in make_folder_list:
    name = name + "_d%s" % (v)
    working_directory = os.path.join(working_dir,
                                     Name_of_project + "_results")
    dest_dir = os.path.join(working_dir,
                            Name_of_project + "_results", name)
    try:
        os.makedirs(dest_dir)
    except OSError:
        print ("folder already exists, I will write over what is in there!!")

###################################################################
# run the program
parse_tab_file_get_clusters(fasta, seq_db, all_fasta,
                            in_file, old_to_new,
                            min_novel_cluster_threshold,
                            right_total_reads,
                            working_directory, v,
                            align, out_file)

