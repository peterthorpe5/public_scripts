#!/usr/bin/env python
#
# tools.py
#
# Reimplementations of Santi's trim_longitudes.py and blastclust_lst2fasta.py
# scripts/functions, and other miscellaneous functions
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os

import gzip
from subprocess import check_output, CalledProcessError
import hashlib
import sys
#import pysam
#import sklearn

from collections import defaultdict
from optparse import OptionParser

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class NotExecutableError(Exception):
    """Exception raised when expected executable is not executable"""
    def __init__(self, message):
        self.message = message


def is_exe(filename):
    """Returns True if path is to an executable file"""
    if os.path.isfile(filename) and os.access(filename, os.X_OK):
        return True
    else:
        try:
            exefile = check_output(["which", filename]).strip()
        except CalledProcessError:
            raise NotExecutableError("{0} does not exist".format(filename))
    return os.path.isfile(exefile) and os.access(exefile, os.X_OK)


# The function below replicates Santi's trim_longitudes.py script, and is
# this functionality is required to ensure compatibility/consistency with
# his old script for testing. Defaults are equivalent to the default settings
# in Santi's script, excepting the filenames, which we now provide directly.
def trim_seq(infname, outfname, lclip=21, rclip=20, minlen=100):
    """Trims all FASTA sequences in infname by lclip and rclip.
    
    lclip and rclip are the 'left'- and 'right'-hand ends, respectively.
    This function writes the results to outfname. Any clipped sequence with
    length < minlen is not written.
    """
    with open(infname, 'r') as fh:
        # Use generators to save memory
        s_trimmed = (s[lclip:-rclip] for s in SeqIO.parse(fh, 'fasta'))
        return SeqIO.write((s for s in s_trimmed if len(s) >= minlen),
                           outfname, 'fasta')

# Commented out PT's changes, which break backwards compatibility!
#def trim_seq(infname, outfname, lclip=53, rclip=0, minlen=100):
#    """Trims all FASTA sequences in infname by lclip and rclip on the
#    'left'- and 'right'-hand ends, respectively and writes the results to
#    outfname. Any clipped sequence with length < minlen is not written.
#
#    Defaults are equivalent to the default settings are for the current ITS1
#    phy primers and the development ITS1 Phy database. PLEASE DONT BLINDLY
#    USE THESE. The filenames, we now provide directly.
#    """
#    # ensure these are int, as they may have been passed as a string.
#    lclip = int(lclip)
#    rclip = int(rclip)
#    with open(infname, 'r') as fh:
#        # rclip needs to be set up to allow for zero values, using logic:
#        # rclip_coordinate = len(s) - rclip
#        # Use generators to save memory
#        s_trim = (s[lclip:(len(s) - rclip)] for s in SeqIO.parse(fh, 'fasta'))
#        return SeqIO.write((s for s in s_trim if len(s) >= minlen),
#                           outfname, 'fasta')


def convert_fq_to_fa(in_file, out_file):
    """ Function to convert fq to fa
    -i the fastq file to be converted
    -o the desired out fasta file name
    the function will auto detect if the file is .gz
    requires biopython
    """
    # if the file is .gz we need to deal with it
    if in_file.endswith(".gz"):
        in_file = gzip.open(in_file, mode='rt', compresslevel=9,
                            encoding=None, errors=None,
                            newline=None)
    SeqIO.convert(in_file, "fastq", out_file, "fasta")


def deduplicate_and_rename(seqlist, vsearch):
    """Removes duplicates from the passed SeqRecords, and replaces sequence
    IDs with the md5 hash of the sequence, suffixed with the number of
    sequences that were replaced. Returns a tuple (sequences, names_old_to_new)
    where names_old_to_new is a list of tuples of (old_name, hash) pairs.

    - sequences   Iterable of Bio.SeqRecord objects
    """
    abundance = defaultdict(int)
    hash_to_seq = defaultdict(str)
    hash_to_name = defaultdict(str)
    new_to_old = defaultdict(list)
    # compile hashed sequence data - the loop is necessary because we need
    # to run through completely once to obtain abundance info. If memory
    # to store this becomes an issue, we might want to try an alternative
    # approach
    for seq in seqlist:
        seqhash = hashlib.md5(str(seq.seq).encode()).hexdigest()
        abundance[seqhash] += 1
        hash_to_seq[seqhash] = seq.seq
        hash_to_name[seqhash] = seq.id
        new_to_old[seqhash].append(seq.id)
    # generate list of sequences, named for the hash, including abundance data,
    # and return with ID:hash lookup
    seqlist = list()
    for name, abundance_val in abundance.items():
        if vsearch:
            seqname = "{0};size={1};".format(hash_to_name[name],
                                             abundance_val)
            seqlist.append(SeqRecord(id=seqname, description="",
                                     seq=hash_to_seq[name]))
        else:
            seqname = "{0}_{1}".format(name, abundance_val)
            seqlist.append(SeqRecord(id=seqname, description="",
                                     seq=hash_to_seq[name]))
    return(seqlist, new_to_old, hash_to_seq)


def dereplicate_name(fasta, database_out, out, vsearch=False):
    """function to dereplicate the seq. Thus generating an
    abundance of that seq. Rename the seq to something
    swarm will work with.
    - fasta     - fasta file of assembled seq
    - database_out   - ouput old name to new file.
    -out        - outfile
    vsearch is for outputting in a specific vsearch reuired format.
    >name;size=6;
    """
    # open files to write to
    fasta_out = open(out, 'w')
    name_out = open(database_out, "w")
    # convert the file to a list of Seqrecord objects
    seqlist = list(SeqIO.parse(fasta, 'fasta'))
    seqlist, new_to_old, hash_to_seq = deduplicate_and_rename(seqlist, vsearch)
    for key, vals in new_to_old.items():
        out_data = "%s\t%s\n" % (key, "\t".join(vals))
        name_out.write(out_data)
    for seq_record in (seqlist):
        SeqIO.write(seq_record, fasta_out, "fasta")
    # close the open files
    fasta_out.close()
    name_out.close()


def check_OTU_db_abundance_val(OTUfasta):
    """this is a function to check that the OTU database has _1
    at the end of the title, this is required for swamr clustering
    takes in a fasta and if required, writes out a new one"""
    # we will iterate through first to see if we need to rewrite it.
    bad_count = 0
    for seq_record in SeqIO.parse(OTUfasta, "fasta"):
        if not seq_record.id.endswith("_1"):
            # print ("OTU db %s missing abundance value" % seq_record.id)
            bad_count = bad_count + 1
    if bad_count > 0:
        # only if we need to write to a new file, we will.
        outfile = OTUfasta.split(".fa")[0] + "abundance.fasta"
        f_out = open(outfile, 'w')
        for seq_record in SeqIO.parse(OTUfasta, "fasta"):
            seq_record.id = seq_record.id + "_1"
            seq_record.description = ""
            SeqIO.write(seq_record, f_out, "fasta")
        f_out.close()
        return outfile
    else:
        return OTUfasta


def filter_sam_file(in_sam, outfile):
    """function for pysam to filter the sam file to obtain matches
    that only match for its length
    # AS:i:0 are perfect matches
    # samfile format, good apge:
    # http://www.metagenomics.wiki/tools/samtools/bam-sam-file-format
    """
    # counter to track number of seq of interest
    no_missmtch = open(outfile, "w")
    samfile = pysam.AlignmentFile(in_sam)
    cig_list = []
    matches = []
    for read in samfile:
        # this apparently are pure matches. But not the same as AS:i:0(?)
        if (len(read.cigar)) == 1:
            out_info = "%s\t%s\t%s\n" % (read.reference_name,
                                         read.query_name,
                                         read.cigar)
            no_missmtch.write(out_info)
            matches.append(out_info)
        cig_list.append(read.cigar)
    no_missmtch.close()
    return cig_list, matches


# Function replacing Santi's blastclust_lst2fasta.py script
def blastclust_to_fasta(infname, seqfname, outdir):
    """Converts input BLASTCLUST output list to a subdirectory of FASTA files.


    Each individual FASTA file contains all sequences from a single cluster.
    The sequences matching the IDs listed in the BLASTCLUST output .lst file 
    should all be found in the same file.

    Returns the output directory and a list of the files, as a tuple.
    """
    outdirname = os.path.join(outdir, "blastclust_OTUs")
    if not os.path.exists(outdirname):
        os.makedirs(outdirname)
    seqdict = SeqIO.index(seqfname, 'fasta')
    outfnames = []
    with open(infname, 'r') as fh:
        otu_id = 0
        for line in fh:
            otu_id += 1
            outfname = os.path.join(outdirname,
                                    "blastclust_OTU_%06d.fasta" % otu_id)
            SeqIO.write((seqdict[key] for key in line.split()),
                        outfname, 'fasta')
            outfnames.append(outfname)
    return (outdirname, outfnames)


# the following four function are to rename the clusters back to their
# original names
def get_names_from_Seq_db(seq_db):
    """function to get a list of name in the seq db"""
    names = []
    names_abudance_removed = []
    for seq_record in SeqIO.parse(seq_db, "fasta"):
        if seq_record.id.endswith("_1"):
            names.append(seq_record.id)
        else:
            names_abudance_removed.append(seq_record.id)
            names.append(seq_record.id + "_1")
    return names, names_abudance_removed


def coded_name_to_species(database_file):
    """functiong takes the already generated tab separated
    database of coded name to species file. Returns a dic
    of coded_name to species"""
    with open(database_file) as file:
        data = file.read().split("\n")
    coded_name_to_species_dict = dict()
    for line in data:
        if not line.strip():
            continue  # if the last line is blank
        if line.startswith("#"):
            continue
        data = line.split("\t")
        coded_name = data[0]
        species = data[1:]
        coded_name_to_species_dict[coded_name.rstrip()] = species
    return coded_name_to_species_dict


def return_real_line(line):
    """function to return only true lines,
    no comments, or blank lines."""
    if not line.strip():
        return False  # if the last line is blank
    if line.startswith("#"):  # dont want comment lines
        return False
    if "\t" in line:
        cluster_line = line.rstrip("\n").split("\t")
    else:
        # different clustering program?
        cluster_line = line.rstrip("\n").split()
    return cluster_line


def parse_tab_file_get_clusters(in_file, seq_db, database, out_file,
                                dev=False):
    """script to open up a tab or space separeted clustering
    output and rename according to the name in the database file.
    Abundance is also appended to the name"""
    # call the function to get the dictionary
    # populated with the database
    names, names_abudance_removed = get_names_from_Seq_db(seq_db)
    coded_name_to_species_dict = coded_name_to_species(database)
    cluster_file = open(in_file, "r")
    summary_out_file = open(out_file, "w")
    count = int(0)
    for line in cluster_file:
        output_str = ""
        count += 1
        # get func to return real line only
        cluster_line = return_real_line(line)
        if not cluster_line:
            continue
        for member in cluster_line:
            if member in names or member in names_abudance_removed:
                if member.endswith("_1"):
                    member = ("_").join(member.split("_")[:-1])
                if dev:  # for development sctripts.PT
                    cluster_summary = "%s_abundance=1\t" % (member)
                else:
                    cluster_summary = "%s\t" % (member)
                output_str = output_str + cluster_summary
                continue
            try:
                # data would be: 2f1454d16278fda2d44f26ebf7a0ed05_310
                # _abundance value
                split_name = member.split("_")[:-1]
                coded_name = ("_").join(split_name)
                if coded_name in names:
                    species = coded_name
                    abundance = "1"
                else:
                    species = coded_name_to_species_dict[coded_name]
                    abundance = member.split("_")[-1]
                    if dev:  # develop for abundance in name
                        species = "\t".join("%s_abundance=%s"
                                            % (i, abundance)
                                            for i in species)
                    else:
                        species = "\t".join("%s" % (i)
                                            for i in species)

            except KeyError:
                print ("we have an error at %s " % member)
                print ("something went wrong w/decoding the names maybe " +
                       "names were not separated in your database? " +
                       "If so, the code above this statement need " +
                       "adjusting accordingly ")
                sys.exit()
            # add the info to a str, we will write at the
            # end of the cluster line
            cluster_summary = "%s\t" % (species)
            output_str = output_str + cluster_summary
        summary_out_file.write(output_str+"\n")
    # close the files
    cluster_file.close()
    summary_out_file.close()


######################################################################
# set of functions to reformat cdhit clusters


def open_parse(clustr):
    """function: opens file and return
    list split on \n from the file"""
    with open(clustr) as file:
        data = file.read().split("\n")
    return data


def get_cluster_member(line):
    """function: split the line from cd hit
    cluster to obtain the member.
    e.g. 0	233nt, >read_1... *
    to return read_1    """
    member = line.split(">")[1]
    member = member.split("...")[0]
    return member


def write_out_clusters(outfile, input_list):
    """funct: write out members of a list to a file.
    The memebers are passed to the function as a list
    [member\t, member_2\t].
    The outfile has already been opened"""
    for i in input_list:
        # should have \t already with the name.
        # but an extra one at the end.
        outfile.write(i)


def reformat_cdhit_clustrs(clustr, outfile, out_R):
    """function: return the cd hit clusters to one line
    per cluster, as is the format from some other clustering
    programs.
    Also, function write the cd hit cluster to a format
    of:
    ID\tcluster_number
    This format is for future R analysis.
    Take in .clstr filen from cdhit
    gets open_parse to return the file as a \n separated list
    write the new format to file."""
    f_out = open(outfile, "w")
    F_R_out = open(out_R, "w")
    data = open_parse(clustr)
    # list to put the cluster memebers in
    cluster_members = []
    R_clusters = []
    for line in data:
        if not line.strip():
            continue  # if the line is blank
        # this is a new cluster group
        if line.startswith(">"):
            # e.g. >Cluster 0
            # starts at zero, need 1 based
            clust_num = int(line.split()[1]) + 1
            cluster_members.append("\n")
            # start of file, memebr list if empty
            if clust_num == 1:
                # remove the first \n
                cluster_members.pop(-1)
                pass
            else:
                pass
        else:
            # these are members of this cluster, can be many!
            # call the func
            member = get_cluster_member(line)
            # this will result in an extra \t after the list
            # entry
            member = member + "\t"
            cluster_members.append(member)
            # format required for R: ID\tcluster_number
            # \n is added to the write out function
            R_format = "%s%d\n" % (member, clust_num)
            R_clusters.append(R_format)
    # call func to write out lists
    write_out_clusters(f_out, cluster_members)
    write_out_clusters(F_R_out, R_clusters)


#################################################
# the following reformat sam for format to work in R


def split_blast_line(line):
    """funct: return a split version of the b;ast6.
    output from vsearch
    Returns: reads, db_entry"""
    elements = line.split("\t")
    reads = elements[0]
    reads = reads.split(";size=")[0]
    db_entry = elements[1]
    db_entry = db_entry.split(";size=")[0]
    return reads.rstrip(), db_entry.rstrip()


def split_sam_line(line):
    """funct: return a split version of the samline.
    Returns: reads, db_entry"""
    if "AS:i:0" not in line:
        print ("WARNING these are not perfect matches")
    elements = line.split("\t")
    reads = elements[0]
    db_entry = elements[2]
    return reads, db_entry


def remove_from_list(remove_id, in_list):
    """funct: remove and item from list.
    But checks it is possible and breaks the program
    if not."""
    try:
        in_list.remove(remove_id)
    except ValueError:
        print ("%s not in passed db file. Exiting!" % remove_id)
        sys.exit(1)
    return in_list


def get_all_seq_names(fasta_in):
    """funct: get all the ids that are used.
    Pass the func the file used that has the db and the reads
    retunrs a list of names"""
    names = []
    for seq_record in SeqIO.parse(fasta_in, "fasta"):
        names.append(seq_record.id)
    return names


def write_out_dict(outfile, outfile2, names, input_dict):
    """funct: write out members of a input_dict to a file.
    The memebers are passed to the function as a dict
    [db_entry] = [read, read2].
    The outfile has already been opened"""
    cluster_count = 0
    for db, reads in sorted(input_dict.items()):
        cluster_count = cluster_count + 1
        data_out = "%s\t%d\n" % (db, cluster_count)
        outfile2.write(db + "\t")
        names = remove_from_list(db, names)
        outfile.write(data_out)
        # now iterate through the reads in the db cluster
        for read in sorted(reads):
            data_out = "%s\t%d\n" % (read, cluster_count)
            outfile.write(data_out)
            outfile2.write(read + "\t")
            names = remove_from_list(read, names)
        outfile2.write("\n")
    # the remaining list should be those which no reads mapped
    # to or singleton reads
    for singleton in sorted(names):
        cluster_count = cluster_count + 1
        data_out = "%s\t%d\n" % (singleton, cluster_count)
        outfile.write(data_out)
        outfile2.write(singleton + "\n")


def reformat_sam_clusters(sam, db_and_reads, outfile):
    """function: return the sam file as
    a dict of [dbname] = read_list
    """
    f_out = open(outfile, "w")
    f_out2 = open(outfile + "_1_line", "w")
    data = open_parse(sam)
    # get list of all db and read names
    names = get_all_seq_names(db_and_reads)
    db_to_reads = defaultdict(list)
    for line in data:
        if not line.strip():
            continue  # if the line is blank
        if line.startswith("@"):
            continue  # these are headers
        reads, db_entry = split_sam_line(line)
        db_to_reads[db_entry].append(reads)
    write_out_dict(f_out, f_out2, names, db_to_reads)
    f_out.close()


def reformat_blast6_clusters(blast6, db_and_reads, outfile):
    """function: return the blast6 file as
    a dict of [dbname] = read_list
    """
    f_out = open(outfile, "w")
    f_out2 = open(outfile + "_1_line", "w")
    data = open_parse(blast6)
    # get list of all db and read names
    names = get_all_seq_names(db_and_reads)
    db_to_reads = defaultdict(list)
    for line in data:
        if not line.strip():
            continue  # if the line is blank
        if line.startswith("	"):
            # this happend if the file is windows formatted
            continue  # these are blast formta
        reads, db_entry = split_blast_line(line)
        db_to_reads[db_entry].append(reads)
    write_out_dict(f_out, f_out2, names, db_to_reads)
    f_out.close()
    f_out2.close()


# the following is to reformat swarm clusters to a format for R
# it calls functions above in this collection
def split_swarm_line(line, dbnames, dbnames_abun_remvd,
                     abundance):
    """funct: return a split version of the swarmline.
    Returns: reads, db_entry.
    Take in a list of db entry names"""
    out_list = []
    elements = line.split()
    for element in elements:
        if element == "":
            continue
        if element.replace("_abundance=", "_") in dbnames:
            # this is a db sequence
            element = element.replace("_abundance=", "_")
        elif element.replace("_abundance=", "_") in dbnames_abun_remvd:
            # this is a db sequence
            element = element.replace("_abundance=", "_1")
        else:
            # this are the assembled reads. They have
            # e.g. _abundance=1 has been added. And need to
            # be removed.
            if abundance:
                element = "".join(element.split("_")[:-1])
            else:
                # this is so we can use blast clust results
                # with the same function
                element = element
        out_list.append(element)
    return out_list


def write_out_swam_to_R(outfile, names, input_lists):
    """funct: write out members of a list of lists.
    The swarm clusters are passed to this function as a list of lists.
    Names is a list of all the sequence- names of which were subjected to
    clustering"""
    cluster_count = 0
    for cluster in input_lists:
        cluster_count = cluster_count + 1
        for member in cluster:
            member = member.rstrip()  # any trailing \n
            data_out = "%s\t%d\n" % (member, cluster_count)
            outfile.write(data_out)
            names = remove_from_list(member, names)


def reformat_swarm_cls(swarm, seq_db, db_and_reads, outfile,
                       abundance=True):
    """function: return the swarm file as
    name\tcluster_number
    This is so we can easily import the data into R.
    Take in a swarm file: 1 line per cluster (\t separated)
    seq_db is the db used to cluster swarm with.
    db_and_reads is a file containing the db and the assmebled
    read."""
    all_clusters = []
    f_out = open(outfile, "w")
    data = open_parse(swarm)
    # get list of all db and read names
    dbnames, dbnames_abun_remvd = get_names_from_Seq_db(seq_db)
    names = get_all_seq_names(db_and_reads)
    for line in data:
        if not line.strip():
            continue  # if the line is blank
        if line.startswith("#"):
            continue  # these are comments
        # This func is needed as the _abundance is not how the other
        # clustering tools will report the names.
        # need to keep the names consitent
        line = line.rstrip()
        cluster = split_swarm_line(line, dbnames, dbnames_abun_remvd,
                                   abundance)
        all_clusters.append(cluster)
    write_out_swam_to_R(f_out, names, all_clusters)
    f_out.close()
