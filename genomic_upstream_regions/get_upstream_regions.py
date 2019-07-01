#!/usr/bin/env python
# TITLE: script to get the upstream regions of genes of interest
# script will return upt to the gene if the full length
# falls within that gene.
# also, script will return reverse complemnet of negative strand coded genes.
# author: Peter Thorpe September 2015.
# The James Hutton Insitute, Dundee, UK.

"""
script to return a user defined threshoold upstream number of nucleotides of
genes of interest from the genome seq using a modified GFF output.
"""

# imports
from Bio.Seq import Seq
from Bio import SeqIO
import time
import os
import errno
import logging
import logging.handlers
import sys
from optparse import OptionParser


def sys_exit(msg, error_level=1):
   """Print error message to stdout and quit with
    given error level."""
   sys.stderr.write("%s\n" % msg)
   sys.exit(error_level)


def wanted_genes(genes_file):
    """function to return a list of wanted genes from file.
    Take in a file which has gen names only. Must match the GFF
    names. Returns a set to avoid duplicates.
    Take in a file.
    Gene names \n sep"""
    wanted = open(genes_file, "r")
    names = wanted.readlines()
    wanted_data = [line.rstrip() for line in names
                   if line.strip() != "" if not line.startswith("#")]
    wanted.close()
    wanted_set = set([])
    for i in wanted_data:
        i = i.replace("id=", "")
        i = i.replace(".t1", "")
        wanted_set.add(i.split(";")[0])
    return wanted_set


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
        error = "coord file fields wrong length, should be 5"
        assert len(gff_info.split("\t")) == 5 , "%s" % error
        gene = gff_info.split("\t")[4]
        if ";" in gene:
            gene = gene.split(";")[0]
        # check each gene only represented once
        assert gene not in gene_list, ("duplicate found %s. Reformat -C file." % gene)
        gene_list.append(gene)
        if gene in coordinate_dict.values():
            print ("repeated line in gff sub file")
            continue
        else:
            scaffold_cordinates = gff_info.split("\t")[:]
            coordinate_dict[gene] = scaffold_cordinates
    return coordinate_dict


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
                                     scaffold])
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


def parse_through_gene_coordinates(coordinate_file):
    """function to return a list of genes and coordinates
    with directions.
    It replaces the ID= from the gene names, as sometime this is
    present and others not."""
    data = open(coordinate_file, "r")
    genes_coordinate = data.readlines()
    genes_coordinate = [line.replace("ID=", "").rstrip() \
                        for line in genes_coordinate
                        if line.strip() != "" \
                        if not line.startswith("#")]
    data.close()
    return genes_coordinate


def down_stream_seq_getter(coordinate_file,
                         genome_sequence,
                         upstream,
                         genes_file,
                         outfile,
                         logger,
                         min_len,
                         user_defined_genic=0):
    """this is a function returns the downstream regions of the list
    of genes of interest - a user defined threshold for the
    number of nucleotides upstream is also used"""
    f= open(outfile, 'w')
    threshold = int(upstream)
    min_len = int(min_len)
    # call the func
    wanted = wanted_genes(genes_file)
    assert len(wanted) == len(set(wanted)), "duplicates in wanted list!"
    logger.info("indexing genome")
    Genome_sequence = SeqIO.index(genome_sequence, "fasta")
    Genome_sequence_time=time.time()
    out = 'import genome file took, %.3f' % (Genome_sequence_time
                                             - start_time)
    logger.info(out)
    coordinates = parse_through_gene_coordinates(coordinate_file)
    # [gene] = scaffold, cordinates
    cordinate_dictionary = index_gene_scaffold_coordinates(coordinate_file)
    # we need to get the gene start and stop information
    # assign the user defined number of genic nucleotides to return to int
    user_defined_genic = int(user_defined_genic)
    for gene_name in wanted:
        coordinate_info = cordinate_dictionary[gene_name]
        # assert len(coordinate_full_info.split("\t"))==5
        contig = coordinate_info[0]
        location_start = coordinate_info[1]
        location_stop = coordinate_info[2]
        direction_of_coding = coordinate_info[3]
        # get the altered start (+) strand, and stop (-) strand.
        # if the user want genic regions, user_defined_genic
        # will bring back X amount of the gene
        DNA_start = int(location_start)
        DNA_stop = int(location_stop)
        error = "are the start stop coordinates are correct?"
        assert DNA_start < DNA_stop, "%s" % error
        downstream = (DNA_stop) + threshold
        negative_strand_downstream = (DNA_start) - threshold
        Genome_seq_record = Genome_sequence[contig]
        length_of_contig = len(Genome_seq_record.seq)
        DNA_region_of_interest = Genome_seq_record.seq
        if downstream > length_of_contig:
            warn = " ".join([gene_name,
                             "downstream region may fall off",
                             "end of contig: ",
                             contig])
            downstream = len(Genome_seq_record.seq)
            logger.warning(warn)
        # test if the "upstream region hits a gene
        # slice up the contig for the region of interest
        # [threshold upstream: genestart]
        # or negative [geneend: genestart plus threshold
        DNA_region_of_interest_downstream_positive \
                    = DNA_region_of_interest[DNA_stop:downstream]
        # for negative strand - we reverse complement it
        DNA_region_of_interest_negative_downstream \
                    = DNA_region_of_interest[(negative_strand_downstream - 1 )
                                             :DNA_start]
        DNA_region_of_interest_negative_downstream \
            = DNA_region_of_interest_negative_downstream.reverse_complement()

        if "-" in direction_of_coding:
            new_start = iterate_coordinate_dict(cordinate_dictionary,
                                                gene_name,
                                                contig,
                                                negative_strand_downstream,
                                                logger,
                                                "downstream")
            if new_start:
                new_start = int(new_start)
                warn = " ".join([gene_name,
                                 "query request is going to return",
                                 "a region of a gene",
                                 "\n",
                                 "going to return upto, New start: ",
                                 str(new_start)])
                logger.warning(warn)

                DNA_region_of_interest_negative_downstream2 = \
                                DNA_region_of_interest[new_start:DNA_start]
                DNA_region_of_interest_negative_downstream = \
                    DNA_region_of_interest_negative_downstream2.reverse_complement()
                if len(DNA_region_of_interest_negative_downstream) > min_len:
                    seq_record = '>%s\t|%s\t[%s:%s]%sbp_downstream - strand\n%s\n' % \
                                 (gene_name,contig,
                                 str(new_start),
                                 DNA_start,
                                 threshold,
                                 DNA_region_of_interest_negative_downstream)
                    f.write(seq_record)
                    if "NNNNN" in DNA_region_of_interest_negative_downstream:
                        warn = "NNN found in %s" % (gene_name)
                        logger.info(warn)
            else:
            # (-) ... downstream is threshold - (DNA_start)
                if len(DNA_region_of_interest_negative_downstream) > min_len:
                    seq_record = '>%s\t|%s\t[%s:%s]%sbp_downstream - strand\n%s\n' % \
                                 (gene_name,
                                  contig,
                                  negative_strand_downstream,
                                  DNA_start,
                                  threshold,
                                  DNA_region_of_interest_negative_downstream)
                    f.write(seq_record)
                    if "NNNNN" in DNA_region_of_interest_negative_downstream:
                        warn = "NNN found in %s" % (gene_name)
                        logger.warning(warn)

        if "+" in direction_of_coding:
            # test if the "upstream region hits a gene
            new_region_to_return = iterate_coordinate_dict(cordinate_dictionary,
                                                           gene_name,
                                                           contig,
                                                           downstream,
                                                           logger,
                                                           "downstream")
            if new_region_to_return:
                new_region_to_return = int(new_region_to_return)
                DNA_region_of_interest_downstream_positive = \
                        DNA_region_of_interest[DNA_stop:int(new_region_to_return)]
                if len(DNA_region_of_interest_downstream_positive) > min_len:
                    seq_record = '>%s\t|%s\t[%s:%s]%sbp_downstream + strand\n%s\n' % \
                                 (gene_name,
                                  contig,\
                                  str(DNA_stop),
                                  str(new_region_to_return),
                                  threshold,
                                  DNA_region_of_interest_downstream_positive)
                    f.write(seq_record)
                    if "NNNNN" in DNA_region_of_interest_downstream_positive:
                        warn = "NNN found in %s" % (gene_name)
                        logger.warning(warn)
            else:
            # (+) ... upstream is (DNA_start) - threshold"
                if len(DNA_region_of_interest_downstream_positive) > min_len:
                    seq_record = '>%s\t|%s\t[%s:%s]%sbp_downstream + strand\n%s\n' % \
                                  (gene_name,
                                  contig,\
                                  str(DNA_stop),
                                  downstream,
                                  threshold,
                                  DNA_region_of_interest_downstream_positive)
                    f.write(seq_record)
                    if "NNNNN" in DNA_region_of_interest_downstream_positive:
                        warn = "NNN found in %s" % (gene_name)
                        logger.warning(warn)
    f.close()


def up_stream_seq_getter(coordinate_file,
                         genome_sequence,
                         upstream,
                         genes_file,
                         outfile,
                         logger,
                         min_len,
                         user_defined_genic=0):
    """this is a function returns the upstream regions of the list
    of genes of interest - a user defined threshold for the
    number of nucleotides upstream is also used"""
    f= open(outfile, 'w')
    threshold = int(upstream)
    min_len = int(min_len)
    # call the func
    wanted = wanted_genes(genes_file)
    assert len(wanted) == len(set(wanted)), "duplicates in wanted list!"
    logger.info("indexing genome")
    Genome_sequence = SeqIO.index(genome_sequence, "fasta")
    Genome_sequence_time=time.time()
    out = 'import genome file took, %.3f' % (Genome_sequence_time
                                             - start_time)
    logger.info(out)
    coordinates = parse_through_gene_coordinates(coordinate_file)
    # [gene] = scaffold, cordinates
    cordinate_dictionary = index_gene_scaffold_coordinates(coordinate_file)
    # we need to get the gene start and stop information
    # assign the user defined number of genic nucleotides to return to int
    user_defined_genic = int(user_defined_genic)
    for gene_name in wanted:
        coordinate_info = cordinate_dictionary[gene_name]
        # assert len(coordinate_full_info.split("\t"))==5
        contig = coordinate_info[0]
        location_start = coordinate_info[1]
        location_stop = coordinate_info[2]
        direction_of_coding = coordinate_info[3]
        # get the altered start (+) strand, and stop (-) strand.
        # if the user want genic regions, user_defined_genic
        # will bring back X amount of the gene
        DNA_start = int(location_start)
        DNA_stop = int(location_stop)
        error = "are the start stop coordinates are correct?"
        assert DNA_start < DNA_stop, "%s" % error
        upstream = (DNA_start) - threshold
        if upstream < 2:
            warn = " ".join([gene_name,
                             "upstream region may fall off",
                             "start of contig: ",
                             contig])
            logger.warning(warn)
        negative_strand_upstream = (DNA_stop) + threshold
        Genome_seq_record = Genome_sequence[contig]
        length_of_contig = len(Genome_seq_record.seq)
        DNA_region_of_interest = Genome_seq_record.seq
        if negative_strand_upstream > length_of_contig:
            warn = " ".join([gene_name,
                             "upstream region may fall off",
                             "end of contig: ",
                             contig])
            logger.warning(warn)
        # test if the "upstream region hits a gene
        # slice up the contig for the region of interest
        # [threshold upstream: genestart]
        # or negative [geneend: genestart plus threshold
        DNA_region_of_interest_upstream_positive \
                    = DNA_region_of_interest[upstream:
                                            (DNA_start +
                                            (user_defined_genic - 1))]
        # for negative strand - we reverse complement it
        DNA_region_of_interest_negative_upstream2 \
                    = DNA_region_of_interest[(DNA_stop - user_defined_genic):
                                            (negative_strand_upstream - 1)]
        DNA_region_of_interest_negative_upstream \
            = DNA_region_of_interest_negative_upstream2.reverse_complement()

        if "-" in direction_of_coding:
            new_start = iterate_coordinate_dict(cordinate_dictionary,
                                                gene_name,
                                                contig,
                                                negative_strand_upstream,
                                                logger,
                                                "upstream")
            if new_start:
                new_start = int(new_start)
                warn = " ".join([gene_name,
                                 "query request is going to return",
                                 "a region of a gene",
                                 "\n",
                                 "going to return upto, New start: ",
                                 str(new_start)])
                logger.warning(warn)

                DNA_region_of_interest_negative_upstream2 = \
                                DNA_region_of_interest[(DNA_stop -
                                                        (user_defined_genic)):
                                                        new_start]
                DNA_region_of_interest_negative_upstream = \
                    DNA_region_of_interest_negative_upstream2.reverse_complement()
                if len(DNA_region_of_interest_negative_upstream) > min_len:
                    seq_record = '>%s\t|%s\t[%s:%s]%sbp_upstream - strand\n%s\n' % \
                                 (gene_name,contig,
                                 (DNA_stop - user_defined_genic),
                                 str(new_start),
                                 threshold,
                                 DNA_region_of_interest_negative_upstream)
                    f.write(seq_record)
                    if "NNNNN" in DNA_region_of_interest_negative_upstream:
                        warn = "NNN found in %s" % (gene_name)
                        logger.warning(warn)
            else:
            # (-) ... downstream is (DNA_stop) + threshold"
                if len(DNA_region_of_interest_negative_upstream) > min_len:
                    seq_record = '>%s\t|%s\t[%s:%s]%sbp_upstream - strand\n%s\n' % \
                                 (gene_name,
                                  contig,
                                 (DNA_stop - user_defined_genic),
                                  negative_strand_upstream,
                                  threshold,
                                  DNA_region_of_interest_negative_upstream)
                    f.write(seq_record)
                    if "NNNNN" in DNA_region_of_interest_negative_upstream:
                        warn = "NNN found in %s" % (gene_name)
                        logger.warning(warn)

        if "+" in direction_of_coding:
            # test if the "upstream region hits a gene
            new_region_to_return = iterate_coordinate_dict(cordinate_dictionary,
                                                           gene_name,
                                                           contig,
                                                           upstream,
                                                           logger,
                                                           "upstream")
            if new_region_to_return:
                new_region_to_return = int(new_region_to_return)
                DNA_region_of_interest_upstream_positive = \
                        DNA_region_of_interest[new_region_to_return:(DNA_start +(
                            user_defined_genic - 1))]
                if len(DNA_region_of_interest_upstream_positive) > min_len:
                    seq_record = '>%s\t|%s\t[%s:%s]%sbp_upstream + strand\n%s\n' % \
                                 (gene_name,
                                  contig,\
                                  str(new_region_to_return),
                                 (DNA_start + user_defined_genic),
                                  threshold,
                                  DNA_region_of_interest_upstream_positive)
                    f.write(seq_record)
                    if "NNNNN" in DNA_region_of_interest_upstream_positive:
                        warn = "NNN found in %s" % (gene_name)
                        logger.warning(warn)
            else:
            # (+) ... upstream is (DNA_start) - threshold"
                if len(DNA_region_of_interest_upstream_positive) > min_len:
                    seq_record = '>%s\t|%s\t[%s:%s]%sbp_upstream + strand\n%s\n' % \
                                  (gene_name,
                                  contig,\
                                  upstream,
                                  (DNA_start + user_defined_genic),
                                   threshold,
                                  DNA_region_of_interest_upstream_positive)
                    f.write(seq_record)
                    if "NNNNN" in DNA_region_of_interest_upstream_positive:
                        warn = "NNN found in %s" % (gene_name)
                        logger.warning(warn)
    f.close()

##########################

if "-v" in sys.argv or "--version" in sys.argv:
    print ("v0.0.4")
    sys.exit(0)


usage = """Use as follows:

$ python get_upstream_regions.py --coordinates
        coordinate_file.fasta -g genome_sequence.fasta
        -upstream <int> number of nucleotides upstream of strat of gene to return
        e.g.  -u 1000
        -z user_defined_genic (how much of the gene to return)
        -m min length of seq to return
        -o outfile_name

Requirements:
python and biopyton.

This will return (--upstream number) of nucleotides to the start of your genes(s) of
interest (-g) gene_file using data from (-c). Gene file can either be space, tab or \n
separated..

The coordinate file can be generated using a GFF3 file and a linux command line using
or using to other python script in this folder. See the readme:

grep "gene" name.gff3 | grep -v "#" | cut -f1,4,5,7,9 > format_for_py_script.out


yeilding this reulting file:

scaffold	start	stop	strand(+/-)	ID=gene_number
GROS_00001	2195	3076	-	ID=GROS_g00002
GROS_00001	8583	10515	+	ID=GROS_g00005.....

The script will check that the start is always less than the end. GFF file should have
starts < stop irrespective of the coding direction

To get all the genes in the file do:

cut -f5 format_for_py_script.out > all_gene_names.out

MORE help / example:

This is an example I ran for G. pallida:

python get_upstream_regions.py -c format_for_py_script.out -g Gpal.v1.0.fas
-f all_gene_names.out -u 225 -z -125
-o Gp.all_gene_names.out_225up_125genic.fasta
> warning_all_gene_names.out_225up_125genic.out


This got (-u 225) 225 bp upstream and (-z 125) 125bp into the current gene for all the
genes in (-f) all_gene_names.out. By default -z is zero. So you dont need to specify
this,
unless you specifically want a piece of the current gene being searched for.
"""

parser = OptionParser(usage=usage)

parser.add_option("-c", "--coordinates",
                  dest="coordinate_file",
                  default="format_for_py_script.out",
                  help="NOTE: coordinate_file can generate using " +
                  "linux command line of "
                  "GFF file:  grep 'gene' name.gff3 | grep -v '#' | " +
                  " cut -f1,4,5,7,9 > format_for_py_script.out ."
                  "Default = format_for_py_script.out")

parser.add_option("-g", "--genome",
                  dest="genome_sequence",
                  default=None,
                  help="genome_sequence.fasta - this has to be the file " +
                  "used to generate the gene models/GFF file")

parser.add_option("-f", "--gene_names",
                  dest="genes_file",
                  default=None,
                  help="a file with a list of gene names to get " +
                  "the upstream regions for")

parser.add_option("-u", "--upstream",
                  dest="upstream",
                  default=False,
                  help="the amount of nucleotide upstream of the gene " +
                  "start, taking into account gene directions, to " +
                  "return in the outfile by default this will not return " +
                  "sequences of min_lenbp or less.")

parser.add_option("-d", "--downstream",
                  dest="downstream",
                  default=False,
                  help="the amount of nucleotide downstream of the gene " +
                  "start, taking into account gene directions, to " +
                  "return in the outfile by default this will not return " +
                  "sequences of min_lenbp or less.")

parser.add_option("-z", "--user_defined_genic",
                  dest="user_defined_genic",
                  default="0",
                  help="the number of nucleotides from within the " +
                  "gene to return, default is 0")

parser.add_option("-m", "--min_len",
                  dest="min_len",
                  default="30",
                  help="the min length of seq to return. " +
                  "Any fragments less than this are not returned " +
                  "Default = 30")

parser.add_option("-o", "--output",
                  dest="out_file",
                  default="upstream_of_genes.fasta",
                  help="Output filename (fasta file)",
                  metavar="FILE")

# get the user options. TODO. arg parser instead
# TODO: Downstream regions.
(options, args) = parser.parse_args()
coordinate_file = options.coordinate_file
genome_sequence = options.genome_sequence
upstream = options.upstream
downstream = options.downstream
genes_file = options.genes_file
outfile = options.out_file
min_len = options.min_len
user_defined_genic = options.user_defined_genic
logfile = outfile.split(".fa")[0] + "WARNINGS.log"

# Run as script
if __name__ == '__main__':
    start_time = time.time()
    # Set up logging
    logger = logging.getLogger('get_upstream_regions.py: %s'
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
    file_list = [coordinate_file, genome_sequence, genes_file]
    for user_file in file_list:
        if not os.path.isfile(user_file):
           print("file not found: %s" % user_file)
           os._exit(0)
    if upstream:
       up_stream_seq_getter(coordinate_file,
                            genome_sequence,
                            upstream,
                            genes_file,
                            outfile,
                            logger,
                            min_len,
                            user_defined_genic)
    if downstream:
       down_stream_seq_getter(coordinate_file,
                              genome_sequence,
                              downstream,
                              genes_file,
                              outfile,
                              logger,
                              min_len,
                              user_defined_genic)
    end_time=time.time()
    logger.info('that took, %.1f' % (end_time - start_time))
