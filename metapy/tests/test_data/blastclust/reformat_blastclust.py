from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys


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


def open_parse(clustr):
    """function: opens file and return
    list split on \n from the file"""
    with open(clustr) as file:
        data = file.read().split("\n")
    return data


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


def write_out_swam_to_R(outfile, names, input_lists):
    """funct: write out members of a list of lists"""
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


reformat_swarm_cls("trimmed.fasta.blastclust99.lst",
                   "phy_db_forblastclust.fasta",
                   "trimmed.fasta",
                   "bc_to_R",
                   False)
