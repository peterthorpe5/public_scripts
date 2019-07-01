#!/usr/bin/env python3
# title: GFF to fasta
# (c) The James Hutton Institute 2018
# Author: Peter Thorpe. The James Hutton Institute, Uk
# imports
import sys
import os
from optparse import OptionParser
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import logging
import logging.handlers
import time

print("""Warning: This is made for a specific purpose. We specific data
      If you are using it, check it works for you!!""")


def index_gene_scaffold_coordinates(coordinate_file):
    """function to return dictionary genes and coordinates
    without directions
    gene = scaffold_cordinates"""
    coordinate_dict = dict()
    data = open(coordinate_file, "r")
    genes_coordinate = data.readlines()
    genes_coordinate = [line.replace("ID=", "").rstrip() \
                        for line in genes_coordinate
                        if line.rstrip() != "" if not line.startswith("#")]
    gene_list = []
    data.close()
    for gff_info in genes_coordinate:
        error = "coord file fields wrong length, should be 9"
        assert len(gff_info.split("\t")) == 9 , "%s" % error
        gene = gff_info.split("\t")[-1]
        gene = gene.replace("ID=gene:", "")
        gene = gene.replace("gene:", "")
        gene = gene.split(";")[0]
        gene = gene.split("Onto")[0]
        if "ncRNA" in gene:
            continue
        # check each gene only represented once
        assert gene not in gene_list, ("duplicate found %s. Reformat -C file." % gene)
        gene_list.append(gene)
        if gene in coordinate_dict.values():
            print ("repeated line in gff sub file")
            continue
        else:
            scaf, source, feature, start, stop, score, \
              direction, frame, gene_info = gff_info.split("\t")
            scaffold_cordinates = [scaf, start, stop, direction, gene]
            coordinate_dict[gene] = scaffold_cordinates
    return coordinate_dict


def index_genome_file(genome):
    """index the genome file
    return a seq_record object that can be accessed
    in a dictionary like manner"""
    genome_database = SeqIO.index(genome, "fasta")
    return genome_database


def check_line(line):
    """checks if line is ok
    Starts with a comment of blank ine = not ok. """
    if line.startswith("#"):
        return False
    if not line.strip():
        return False
    return line


def split_line(line):
    """split the gff line
    Takes in a gff line. returns the elements, as str or int"""
    warning_list = ["gene", "exon", "intron"]
    assert len(line.split("\t")) == 9 , "GFF fields wrong length should be 9"
    scaff, source, feature, start, stop, score, \
            direction, frame, gene_info = line.split("\t")
    if feature in warning_list:
        print("This script is for transtart output only. Be carfeful!!!")
    gene_info = gene_info.rstrip()
    start = int(start)
    stop = int(stop)
    return scaff, source, feature, start, stop, score, \
            direction, frame, gene_info



def iterate_coordinate_dict(coordinate_dict,
                            gene_name,
                            scaffold,
                            coordinates,
                            logger,
                            direction):
    """check to see if the scaffold and new coordinate hits a
    predicted gene. It wanrs the user if the desired upstream
    regions falls into an existing gene. The coordinates will
    be altered upto this gene that is hits."""
    for gene, vals in coordinate_dict.items():
        # find the genes on the same scaffold
        if scaffold in vals[0]:
            dictionary_scaffold = vals[0]
            gene = vals[4]
            # if its is the same gene as the stop
            if gene_name == gene:
                continue
            if scaffold != dictionary_scaffold:
                continue
            else:
                # debugging comment due to Roman numeral scaffold
                # name being "within" eachother
                # print ("scaffold = ", scaffold,
                #        "acutally looking at", vals[0])
                start = int(vals[1])
                stop = int(vals[2])
                coordinates = int(coordinates)
                direction = vals[3]
                # print ("coord=%s, start=%s, stop=%s, gene_query=%s,
                        # gene_GFF=%s" %(coordinates,
                        # start, stop, gene_name, gene))
                # basically does the coordinate fall in the current
                # coordinate for a gene
                if coordinates >= start and coordinates <= stop:
                    warn = " ".join([gene_name,
                                     "upstream coordinate falls in the",
                                     "genic regions of",
                                     gene,
                                     "on scaffold",
                                     scaffold,
                                     str(vals)])
                    logger.warning(warn)
                    data = coordinate_dict[gene_name]
                    if direction == "upstream":
                       if "+" in data:
                           # + coding gene, upstream will be
                           # returned as the end (stop)
                           # of the preceding gene
                           return stop
                       else:
                           # - coding gene, upstream will be
                           # returned as the begining (start)
                           # (stop) of the proceding gene
                           return start
                    if direction == "downstream":
                       if "+" in data:
                           # + coding gene, downstream will be
                           # returned as the star
                           # of the proceding gene
                           return start
                       else:
                           # - coding gene, downstream will be
                           # returned as the begining (stop)
                           # (stop) of the preceding gene
                           return stop
    return False



def gff_to_fasta(gff, genome, coordinate_dict, min_length, Max_length,
                 outfile, upstream, into_TSS, NNN):
    """take in gff file. Gets the seq defined by the gff coords.
    If negative direction coding, the reverse complement is generated.
    A min length of seq to return and max len is applied to remove seq
    less than, for example 3 which cant be real and less that e.g., 25k
    which will be flase positives and not informative in downstream analysis
    """
    print("Indexing the genome")
    min_length = int(min_length)
    genome_database = index_genome_file(genome)
    print("Now iterating through the GFF. Assume it is sorted")
    f_out = open(outfile, "w")
    bind_out = outfile.split(".fa")[0] + "_%dnt_upstream_%d_into_TSS.fasta" % (upstream,
                                                                               into_TSS)
    bind_out_fa = open(bind_out, "w")
    upstream = int(upstream)
    NNN_reject_count = 0
    starting_UTR_count = 0
    fa_out_count = 0
    missing = 0
    with open(gff, "r") as f_handle:
        for line in f_handle:
            line = check_line(line)
            if not line:
                continue
            scaff, source, feature, start, stop, score, \
            direction, frame, gene_info = split_line(line)
            gene_info = gene_info.split("Onto")[0]
            seq_record = genome_database[scaff]
            starting_UTR_count += 1
            if direction == "+":
                bind_seq = seq_record.seq[(start - upstream):(start + into_TSS)]
                out = "\t".join([gene_info, "start ;", str(start), direction, "adjust new values: ", str(upstream), str(into_TSS)])
                logger.info(out)
                new_co_stop = start - upstream
                new_co_start = start + into_TSS
                out = "\t".join(["new_co_start ;", str(new_co_start), "new_co_stop ;", str(new_co_stop), "+"])
                logger.info(out)

                UTR = seq_record.seq[new_co_start:new_co_stop]
                # check the new coordinate does not hit a gene
                new_start = iterate_coordinate_dict(coordinate_dict,
                                                    gene_info,
                                                    scaff,
                                                    new_co_start,
                                                    logger,
                                                    "upstream")
                out = "\t".join(["new_start ; ", str(new_start),  direction ])
                logger.info(out)

                if new_start:
                    new_start = int(new_start)
                    warn = " ".join([gene_name,
                                     "query request is going to return",
                                     "a region of a gene",
                                     "\n",
                                     "going to return upto, New start: ",
                                     str(new_start)])
                    logger.warning(warn)
                    # assign a new region for the seq record . seq
                    UTR = seq_record.seq[new_start:new_co_stop]
                    out = "\t".join(["new_start ;", str(new_start), "new_co_stop ;", str(new_co_stop)])
                    logger.info(out)
            if direction == "-":
                UTR = reverse_complement(seq_record.seq[start:stop])
                bind_seq = reverse_complement(seq_record.seq[(stop - into_TSS)
                                                             :(stop + upstream)])
                new_co_start = stop - upstream
                new_co_stop = stop + into_TSS
                UTR = reverse_complement(seq_record.seq[new_co_start:new_co_stop])
                out = (gene_info, "start ;", start, direction, "adjust new values: ", upstream, into_TSS)
                logger.info(out)
                out = "\t".join(["new_co_start ;", str(new_co_start), "new_co_stop ;", str(new_co_stop), "-"])
                logger.info(out)
                ##################################
                new_stop = iterate_coordinate_dict(coordinate_dict,
                                                    gene_info,
                                                    scaff,
                                                    new_co_stop,
                                                    logger,
                                                    "upstream")
                if new_stop:
                    new_stop = int(new_stop)
                    warn = " ".join([gene_name,
                                     "query request is going to return",
                                     "a region of a gene",
                                     "\n",
                                     "going to return upto, new_stop: ",
                                     str(new_stop)])
                    logger.warning(warn)
                    # assign a new region for the seq record . seq
                    UTR = reverse_complement(seq_record.seq[new_co_start:new_stop])
                    UTR = seq_record.seq[new_start:new_stop]
                    logger.info("new_start ;", new_start, "new_stop ;", new_stop)
                ###########
            coordinate_info = "\tScaffold: %s UTR_and_TSS: %d:%d  Coding_direction: %s  " % (scaff,
                                                                                 start,
                                                                                 stop,
                                                                                 direction)
            fastaextra = "returning: %d upstream of TSS and %d into UTR and or gene:" % (upstream,
                                                                                     into_TSS)
            new_coord = "  %d: %d" % (new_co_start, new_co_stop)
            description = coordinate_info + fastaextra + new_coord
            outstr = ">%s\n%s\n" % (gene_info, UTR)
            bind_str = ">%s_TSS%s\n%s\n" % (gene_info,
                                            description,
                                            bind_seq)
            if len(UTR) > min_length and len(UTR) < Max_length:
                f_out.write(outstr)
                if (NNN*"N") in bind_seq:
                    NNN_reject_count += 1
                    continue  #  we dont want NNNs
                if len(bind_seq) >= upstream: 
                    bind_out_fa.write(bind_str)
                    fa_out_count += 1
            else:
                missing += 1
    print("%d number of genes were rejected due to NNN in seq" % NNN_reject_count)
    print("we had %d UTR and TSS predictions" % starting_UTR_count)
    print("we have outputted %d fasta sequences" % fa_out_count)
    print("missing due to length problems %d" % missing)
    f_out.close()
    bind_out_fa.close()


#############################################################################
#to run it:


usage = """Use as follows:

python GFF_to_fasts.py --gff transtart.gff -m min len -x max len
        -g genome.fasta -o UTR.fasta

"""

parser = OptionParser(usage=usage)

parser.add_option("--gff", dest="gff",
                  default=None,
                  help="the predicted coordinate for the cds predictions .gff",
                  metavar="FILE")

parser.add_option("--genome_gff", dest="genome_gff",
                  default=None,
                  help="the predicted genes as genome_gff.gff",
                  metavar="FILE")

parser.add_option("-g", "--genome",
                  dest="genome",
                  default=None,
                  help="the genome sequence.",
                  metavar="FILE")

parser.add_option("-m", "--min_length",
                  dest="min_length",
                  default=8,
                  help="min_length of seq to return",
                  metavar="FILE")

parser.add_option("-x", "--max_length",
                  dest="max_length",
                  default=4000,
                  help="max_length of seq to return",
                  metavar="FILE")

parser.add_option("-i", "--into_TSS",
                  dest="into_TSS",
                  default=0,
                  help="into_TSS of seq to return to the binding theory",
                  metavar="FILE")

parser.add_option("-u", "--upstream",
                  dest="upstream",
                  default=20,
                  help="upstream of TSS to return",
                  metavar="FILE")

parser.add_option("-n", "--NNN",
                  dest="NNN",
                  default=300,
                  help="number of NNN in a row to be rejected from the output",
                  metavar="FILE")

parser.add_option("-o", "--out", dest="outfile",
                  default="transtart.UTR.fasta",
                  help="Output filename (transtart.UTR.fasta)",
                  metavar="FILE")

(options, args) = parser.parse_args()

#-g
genome = options.genome
#--gff
gff = options.gff
#-o
outfile= options.outfile
# - upstream
upstream = int(options.upstream)
# -x
max_length = int(options.max_length)
# --min_length
min_length = int(options.min_length)
# into_TSS
into_TSS = int(options.into_TSS)
# -n
NNN = int(options.NNN)
logfile = outfile.split(".fa")[0] + "WARNINGS.log"

#######################################################################
# Run as script
if __name__ == '__main__':
    # no logging for this.
    if not os.path.isfile(genome):
        sys.exit("Input genome file not found: %s" % genome)
    if not os.path.isfile(gff):
        sys.exit("Input gff file not found: %s" % gff)
    if not os.path.isfile(options.genome_gff):
        sys.exit("Input gff file not found")
    # Run as script
    start_time = time.time()
    # Set up logging
    logger = logging.getLogger('GFF_to_fasta: %s'
                               % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)
    try:
        logstream = open(logfile, 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        logger.error("Could not open %s for logging" %
                     logfile)
        sys.exit(1)
    coordinate_dict = index_gene_scaffold_coordinates(options.genome_gff)
    gff_to_fasta(gff, genome, coordinate_dict, min_length,
                 max_length, outfile, upstream,
                 into_TSS, NNN)
