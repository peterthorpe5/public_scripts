#!/usr/bin/env python
#author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

#imports
import os
import sys
from sys import stdin,argv
import sys
from optparse import OptionParser
from coverage import *
import numpy
import subprocess
import datetime

#Title:
#script to open gff and create a dictionary of {scaffold: set([gene1, gene2])"
# this can then be used to see if all genes on a scaff are predicted to be HGT and therefor
# the whole scaffold is contamination. 

#############################################################################
#functions
##    try:
##        diamond_tab_as_list = read_diamond_tab_file(diamond_tab_output)
##    except IOError as ex:
##        print("sorry, couldn't open the file: " + ex.strerror + "\n")
##        print ("current working directory is :", os.getcwd() + "\n")
##        print ("files are :", [f for f in os.listdir('.')])


#make a temp_folder_for_all_the_out_files
file_name = 'test.txt'
working_dir = os.getcwd()
dest_dir = os.path.join(working_dir, 'temp')
try:
    os.makedirs(dest_dir)
except OSError:
    print ("folder already exists, I will write over what is in there.")

#############################################################################

def split_gff_gene_names(gene_info):
    "function to hopefully return just the gene name"
    #print gene_info
    try:
        gene = gene_info.split("ID=")[1]
    except:
        gene = gene_info.split("Parent=")[1]
    gene = gene.split(".gene")[0]
    gene = gene.split(".exon")[0]
    gene = gene.split(";")[0]
    gene = gene.split(";")[0]
    gene = gene.split(".CDS")[0]
    gene = gene.rstrip("\n")
    return gene

def parse_gff(gff):
    "function to parse GFF and produce a scaffold_to_gene_dict"
    f_in = open(gff, "r")
    # data to output of function
    scaffold_to_gene_dict = dict()
    gene_to_exon_count = dict()
    gene_start_stop_dict = dict()
    count = 1
    #iterate through gff
    for line in f_in:
        if line.startswith("#"):
            continue
        scaffold,aug,cds_type,start,stop,e,f,g,gene_info = line.split("\t")
        gene = split_gff_gene_names(gene_info)
        #print gene

        #scaffold_to_gene_dict
        if line.split("\t")[2] == "gene":
            if not scaffold in scaffold_to_gene_dict:
                scaffold_to_gene_dict[scaffold]=[gene]
            else:
                scaffold_to_gene_dict[scaffold].append(gene)
            start_stop_formatted = "%s\t%d\t%d" %(scaffold, int(start), int(stop))
            gene_start_stop_dict[gene] = start_stop_formatted

        #gene_to_exon_count
        if line.split("\t")[2] == "exon":
            if not gene in gene_to_exon_count:
                gene_to_exon_count[gene] = count
            else:
                gene_to_exon_count[gene]= count +1
    #print scaffold_to_gene_dict
    f_in.close()
    return scaffold_to_gene_dict, gene_to_exon_count, gene_start_stop_dict


def LTG_file(LTG):
    """function to parse LTG prediction
    get a set of names, and a gene_to_comment_dict"""
    # in file is the output of the LTG prediction tool
    f_in = open(LTG, "r")
    # is the gene >70% identical to it's blast hit?
    # if so, maybe a contamination?
    gene_to_comment_dict = dict()
    HGT_predicted_gene_set = set([])
    for line in f_in:
        if line.startswith("#"):
            continue
        gene = line.split("\t")[0]
        comment = line.split("\t")[-1]
        kingdom = line.split("\t")[5]
        gene = gene.replace("ID=", "").split("gene=g")[0]
        gene = gene.rstrip()
        if not ".t" in gene:
            gene = gene.replace("t", ".t")
        HGT_predicted_gene_set.add(gene)
        data_out_formatted = "%s\t%s" %(comment.rstrip("\n"), kingdom)
        gene_to_comment_dict[gene] = data_out_formatted
    #print HGT_predicted_gene_set
    f_in.close()
    return HGT_predicted_gene_set, gene_to_comment_dict

def get_stats_on_AT_content(dna_file):
    """function to get the mean at standard dev for AT content across
    all genes"""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    import numpy
    gene_AT_cont_dic = dict()

    AT_content_list = []

    for seq_record in SeqIO.parse(dna_file, "fasta"):
        seq_record.seq = seq_record.upper()
        sequence = str(seq_record.seq)
        a_count = sequence.count('A')
        t_count = sequence.count('T')
        #count AT content
        AT = a_count+t_count/len(seq_record.seq)
        #put that in a list 
        AT_content_list.append(AT)
        # assign AT to the gene name for testing later
        gene_AT_cont_dic[seq_record.id] = AT
    #calc average AT for all gene
    the_mean = sum(AT_content_list) / float(len(AT_content_list))
    #calc SD for AT for all genes
    standard_dev = numpy.std(AT_content_list)
    return gene_AT_cont_dic, the_mean, standard_dev

def parse_rnaseq(rnaseq):
    """parse the rnaseq-file
    Take data in like this:
    # Name	Length	TPM	NumReads
    g1.t1	3333	82.2211	186993. May need to be modified to suit
    program output of choice. """

    gene_to_expression = dict()
    with open(rnaseq, "r") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            Name, Length, TPM, NumReads = line.rstrip("\n").split()
            data_formatted = "%s\tNumReads = %s" %(TPM, NumReads)
            gene_to_expression[Name]= data_formatted
    return gene_to_expression


def check_HGT_AT_vs_global_AT(gene_AT_cont_dic, the_mean, standard_dev,
                              gene_of_interest, comment, sd_numbers, gene_to_expression,
                              gene_to_exon_count):
    """function to check the AT content of a gene of interest vs the global
    AT using the mean and SD already generated"""
    # user defined number of standard deviations away from the mean for the stats
    sd_numbers = float(sd_numbers)
    current_gene_AT = gene_AT_cont_dic[gene_of_interest]
    #print "gene_of_interest %s has AT cont of %d" %(gene_of_interest, current_gene_AT)
    lower_threshold = float(the_mean) - (sd_numbers*float(standard_dev))
    upper_threshold = float(the_mean) + (sd_numbers*float(standard_dev))
    #print lower_threshold, upper_threshold
    
    if current_gene_AT < lower_threshold or current_gene_AT > upper_threshold:
        #call dict to get expression
        TPM = gene_to_expression[gene_of_interest]
        exons = gene_to_exon_count[gene_of_interest]
        print "gene: %s\tAT_cont: %d\tcomment: ...%s... \texpression: %s\texons: %d"\
                %(gene_of_interest, current_gene_AT, comment, TPM, exons)
        

# main function

def check_scaffolds_for_only_HGT_genes(genome, gff, LTG, dna_file, sd_numbers, rnaseq,
                                       bam_file, out_file):
    """main function. This calls the other function to get a dictionary
    of genes on scaffolds. A list of HGT/LTG genes and check the scaffolds
    to identify those that only have HGT genes on them. If so, then this
    is most likely a contaminant contig/scaffold"""
    out = open(out_file, "w")
    
    #call function to get the scaffold to gene dict
    scaffold_to_gene_dict, gene_to_exon_count, gene_start_stop_dict = parse_gff(gff)
           
    #call function to get gene_set, gene_to_comment_dict
    HGT_predicted_gene_set, gene_to_comment_dict = LTG_file(LTG)

    #get_scaffold_coverage from import coverage.py
    if bam_file:
        overall_coverage = "overall_coverage.txt"
        overall_expression_dic = get_total_coverage(bam_file, overall_coverage)
        scaffold_coverage_dict, mean_scaf_cov, sd_scaf_cov = get_scaffold_coverage(genome, \
                                                            gene_start_stop_dict,\
                                                            bam_file, overall_expression_dic,\
                                                            HGT_predicted_gene_set)

    #print gene_to_exon_count

    #call function to get rna seq mapping TPM
    gene_to_expression = parse_rnaseq(rnaseq)

    #call function with DNA file
    gene_AT_cont_dic, the_mean, standard_dev = get_stats_on_AT_content(dna_file)
    print "the AVR AT = %f with SD %f " %(the_mean, standard_dev)
    
    for gene, comment in gene_to_comment_dict.items():
        check_HGT_AT_vs_global_AT(gene_AT_cont_dic, the_mean,
                                  standard_dev, gene, comment,
                                  sd_numbers, gene_to_expression, gene_to_exon_count)
        

    for scaffold, genes in scaffold_to_gene_dict.items():
        bad_contig = True
        for gene in genes:
            #print gene
            if gene not in HGT_predicted_gene_set:
                bad_contig = False

        if bad_contig == True:
            data_formatted = "%s\tBad scaffold\n" %(scaffold)
            print "Bad scaffold = %s" %(scaffold)
            out.write(data_formatted)
    out.close()


#################################################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

$ python ~/misc_python/Lateral_gene_transfer_prediction_tool/check_contaminants_on_contigs.py --gff ../augustus.gff3 -LTG LTG_LGT_candifates.out (default)

Requires:
samtools 1.2 or later for Bam file - default is not to try this
Biopython
NUmpy


You may have to tidy and sort your GFF to a GFF3. Use GenomeTools
1) convert augustus.gft to gff3
gt gtf_to_gff3 -o test.gff3 -tidy augustus.gtf
or
gt gff3 -sort -tidy augustus.gff > formatted.gff3


"""

parser = OptionParser(usage=usage)

parser.add_option("--gff", dest="gff", default="test.gff",
                  help="hintsfile",
                  metavar="FILE")
parser.add_option("--LTG", dest="LTG", default="test_LTG_LGT_candifates.out",
                  help="LTG outfile. ",
                  metavar="FILE")
parser.add_option("--dna", dest="dna", default=None,
                  help="predicted cds nucleotide genes for AT content stats ",
                  metavar="FILE")
parser.add_option("-g", "--genome", dest="genome", default=None,
                  help="genome.fasta ",
                  metavar="FILE")
parser.add_option("-s", dest="sd_numbers", default=3,
                  help="the number of stadard deviations away from the mean"
                  " for identifying genes "
                  " that differ from normal AT content. default=3")
parser.add_option("-r", "--rnaseq", dest="rnaseq", default=None,
                  help="RNAseq expression profile for genes. "
                  " in format # Name	Length	TPM	NumReads ")
parser.add_option("-b", "--bam", dest="bam_file", default=False,
                  help="RNAseq expression profile for genes. "
                  "bam file (sorted) with genomic reads mapped to geneome "
                  " this is used to see if HGT genes have a different "
                  " genomic coverage compared to other gene. Requires "
                  " samtools 1.2 or later ")

parser.add_option("-o", "--out_file", dest="out_file", default="test_output.txt",
                  help="outfile to list the bad contigs")


(options, args) = parser.parse_args()


gff = options.gff
out_file = options.out_file
LTG = options.LTG
dna = options.dna
sd_numbers = options.sd_numbers
rnaseq = options.rnaseq
bam_file = options.bam_file
genome = options.genome


#run the program

if not os.path.isfile(gff):
    sys_exit("Input gff file not found: %s" % gff)

if not os.path.isfile(LTG):
    sys_exit("Input LTG file not found: %s" % LTG)
if not os.path.isfile(dna):
    sys_exit("Input dna file not found: %s" % dna)
    
    
check_scaffolds_for_only_HGT_genes(genome, gff, LTG, dna, sd_numbers, rnaseq, bam_file, out_file)


