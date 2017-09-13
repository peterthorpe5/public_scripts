#!/usr/bin/env python
# author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

# imports
import os
import sys
from sys import stdin,argv
import sys
from optparse import OptionParser
from coverage import *
import numpy
import subprocess
import datetime

import datetime
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy

# Title:
# script to open gff and create a dictionary of {scaffold: set([gene1, gene2])"
# this can then be used to see if all genes on a scaff are predicted to be
# HGT and therefor
# the whole scaffold is contamination. 

#############################################################################

# make a temp_folder_for_all_the_out_files
file_name = 'test.txt'
working_dir = os.getcwd()
dest_dir = os.path.join(working_dir, 'temp')
try:
    os.makedirs(dest_dir)
except OSError:
    print ("folder already exists, I will write over what is in there.")

#############################################################################

def split_gff_gene_names(gene_info):
    """function to return just the gene name"""
    #print gene_info
    try:
        gene = gene_info.split("ID=")[1]
    except:
        gene = gene_info
    try:
        gene = gene_info.split("Parent=")[1]
    except:
        gene = gene
    gene = gene.split(".gene")[0]
    gene = gene.split(".exon")[0]
    gene = gene.split(";")[0]
    gene = gene.split(";")[0]
    gene = gene.split(".CDS")[0]
    gene = gene.split(".t")[0]
    # some data set require this to be uncommented
    # gene = gene.split(".t")[0]
    gene = gene.rstrip("\n")
    #gene = gene.replace(".t1", "")
    return gene

def parse_gff(gff):
    """function to parse GFF and produce a scaffold_to_gene_dict
    3 dictionaries are returned:
    scaffold_to_gene_dict, gene_to_exon_count, gene_start_stop_dict"""
    f_in = open(gff, "r")
    # data to output of function
    scaffold_to_gene_dict = dict()
    gene_to_exon_count = dict()
    gene_start_stop_dict = dict()
    count = 1
    # iterate through gff
    for line in f_in:
        if line.startswith("#"):
            continue
        if not line.strip():
                continue  #  if the last line is blank
        scaffold, aug, cds_type, start, stop, e, f, \
                  g,gene_info = line.split("\t")
        gene = split_gff_gene_names(gene_info)
        # print gene

        # scaffold_to_gene_dict
        if line.split("\t")[2] == "gene":
            if not scaffold in scaffold_to_gene_dict:
                scaffold_to_gene_dict[scaffold]=[gene]
            else:
                scaffold_to_gene_dict[scaffold].append(gene)
            start_stop_formatted = "%s\t%d\t%d" %(scaffold, int(start),
                                                  int(stop))
            gene_start_stop_dict[gene] = start_stop_formatted

        # gene_to_exon_count
        if line.split("\t")[2] == "exon":
            if not gene in gene_to_exon_count:
                gene_to_exon_count[gene] = 1
            else:
                gene_to_exon_count[gene] += 1
    # print scaffold_to_gene_dict
    f_in.close()
    print("if this crashes check exons are in your gff")
    print("if not change lines 96 'exon' to 'intron'")
    return scaffold_to_gene_dict, gene_to_exon_count,\
           gene_start_stop_dict


def LTG_file(LTG):
    """function to parse LTG prediction
    get a set of names, and a gene_to_comment_dict"""
    # in file is the output of the LTG prediction tool
    f_in = open(LTG, "r")
    # is the gene >70% identical to it's blast hit?
    # if so, maybe a contamination?
    gene_to_HGT_percent_identity = dict()
    gene_to_comment_dict = dict()
    HGT_predicted_gene_set = set([])
    gene_to_HGTspeces_discription_dict = dict()
    gene_to_AI = dict()
    for line in f_in:
        if line.startswith("#"):
            continue
        if not line.strip():
            continue  #  if the last line is blank
        gene = line.split("\t")[0]
        # call function to format gene name
        gene = split_gff_gene_names(gene)
        HGT_percent_identity = line.split("\t")[1]
        comment = line.split("\t")[-1]
        species = line.split("\t")[-3]
        description = line.split("\t")[-2]
        AI = line.split("\t")[7]
        kingdom = line.split("\t")[5]
        gene = gene.replace("ID=", "").split("gene=g")[0]
        gene = gene.rstrip()
        if not ".t" in gene:
            gene = gene.replace("t", ".t")
        HGT_predicted_gene_set.add(gene)
        data_out_formatted2 = "%s\t%s" %(species, description)
        data_out_formatted = "%s\t%s" %(comment.rstrip("\n"),kingdom)
        gene_to_HGTspeces_discription_dict[gene] = data_out_formatted2
        gene_to_comment_dict[gene] = data_out_formatted
        gene_to_AI[gene] = AI
        gene_to_HGT_percent_identity[gene] = HGT_percent_identity
    # print (HGT_predicted_gene_set)
    f_in.close()
    return HGT_predicted_gene_set, gene_to_comment_dict,\
           gene_to_HGTspeces_discription_dict, gene_to_AI, \
           gene_to_HGT_percent_identity


def get_stats_on_AT_content(dna_file):
    """function to get the mean at standard dev for AT content across
    all genes"""

    gene_AT_cont_dic = dict()
    AT_content_list = []
    for seq_record in SeqIO.parse(dna_file, "fasta"):
        seq_record.seq = seq_record.upper()
        sequence = str(seq_record.seq)
        a_count = sequence.count('A')
        t_count = sequence.count('T')
        ##count AT content
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
    program output of choice.
    returns a dictionary of gen to expression values"""

    gene_to_expression = dict()
    with open(rnaseq, "r") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            #print line
            Name, Length, TPM, NumReads = line.rstrip("\n").split()
            #Name = Name.replace(".t1", "")
            ##call function to format gene name
            #Name = split_gff_gene_names(Name)
            data_formatted = "%s\tNumReads = %s" %(TPM, NumReads)
            gene_to_expression[Name]= data_formatted
    return gene_to_expression


def how_many_sd_from_mean(mean, sd, AT_content):
    "function to retunr the number of sd a value is from the mean"
    return (float(mean) - AT_content)/float(sd)
    

def check_genomic_cov(mean, sd, genomic_cov):
    "function to check the genomic cov of gene of interest"

    NUmber_of_sd_from_mean = how_many_sd_from_mean(mean, sd, genomic_cov)
    return NUmber_of_sd_from_mean
    

def check_HGT_AT_vs_global_AT(gene_AT_cont_dic, AI, the_mean, standard_dev,
                              gene_of_interest, comment, sd_numbers,
                              gene_to_expression,
                              gene_to_exon_count,
                              gene_to_HGTspeces_discription_dict,
                              genomic_cov_from_mean,
                              gene_to_HGT_percent_identity):
    """function to check the AT content of a gene of interest vs the global
    AT using the mean and SD already generated"""
    
    # user defined number of standard deviations away from
    # the mean for the stats
    sd_numbers = float(sd_numbers)
    # print "gene_of_interest = ", gene_of_interest
    current_gene_AT = gene_AT_cont_dic[gene_of_interest]

    AI = float(AI)

    num_sd_from_mean = how_many_sd_from_mean(the_mean, standard_dev,
                                             current_gene_AT)
    assert  how_many_sd_from_mean(10, 2, 2) ==4
    HGTspecies_description = gene_to_HGTspeces_discription_dict[gene_of_interest]
    HGTspecies = HGTspecies_description.split("\t")[0]
    description = HGTspecies_description.split("\t")[1]
    # description = description.split("[")[0]
    lower_threshold = float(the_mean) - (sd_numbers * float(standard_dev))
    upper_threshold = float(the_mean) + (sd_numbers * float(standard_dev))
    # print lower_threshold, upper_threshold
    if current_gene_AT < lower_threshold or current_gene_AT > upper_threshold:
        # if calling with RNAseq assembly
        transcript = gene_of_interest.split("|")[0]
        # call dict to get expression
        try:
            TPM = gene_to_expression[gene_of_interest]
        except:
            ValueError #Rnaseq assebmly not genome
            TPM = gene_to_expression[transcript]
        try:
            exons = gene_to_exon_count[gene_of_interest]
        except:
            ValueError
            # note if this crashes you may need to specifcy introns,
            # not exons to count.
            # exons = gene_to_exon_count[transcript]
            exons = 1
        try:
            per_ident = gene_to_HGT_percent_identity[gene_of_interest]
        except:
            ValueError
            per_ident = gene_to_HGT_percent_identity[transcript]
        #print "gene: %s\tAT_cont: %d\tcomment: ...%s... \texpression: %s\texons: %d\t%s\t%s"\
                #%(gene_of_interest, current_gene_AT, comment, TPM, exons, HGTspecies, description)
        print ("gene: %s\tAI = %s\tAT_cont: %d\tAT_cont_numSD_fromMean: %.2f\tgenomic_cov_from_mean: %.2f\t\texpression: %s\texons: %d\tcomment: ...%s...\t%s\t%s" %(gene_of_interest,
                      AI,  current_gene_AT,\
                      num_sd_from_mean, genomic_cov_from_mean, \
                      TPM, exons, comment, HGTspecies, description))
        HGT_info_dataformatted = "%s\t%s\t%.1f\t%d\t%.2f\t%.2f\t%s\t%d\t%s\t%s\t%s\n" %(gene_of_interest,
                                                                                        per_ident,
                                                                                        AI,
                                                                                        current_gene_AT,
                                                                                        num_sd_from_mean,
                                                                                        genomic_cov_from_mean,
                                                                                        TPM,
                                                                                        exons,
                                                                                        comment,
                                                                                        HGTspecies,
                                                                                        description)
    return HGT_info_dataformatted

        

# main function

def check_scaffolds_for_only_HGT_genes(genome, gff, LTG, dna_file, sd_numbers, rnaseq,
                                       bam_file, out_file):
    """main function. This calls the other function to get a dictionary
    of genes on scaffolds. A list of HGT/LTG genes and check the scaffolds
    to identify those that only have HGT genes on them. If so, then this
    is most likely a contaminant contig/scaffold"""
    bad_scaffold_out = "bad_scaffold."+out_file
    bad_scaff_title = "#%s\n#scaffold\tcomment\tGenes_on_scaffold\tAI_of_these_gene\tdescription\n" %(datetime.date.today())
    out = open(bad_scaffold_out, "w")
    out.write(bad_scaff_title)

    HGT_gene_info = "HGT.info."+out_file
    f_out = open(HGT_gene_info, "w")
    f_out.write("#%s\n#gene\tPerc_identity_to_HGT_hit\tAI\tAT_cont\tAT_cont_numSD_fromMean\tgenomic_cov_from_mean\texpression_TMM\tnum_RNAseq_reads\texons\tcomment\tKingdom\tHGT_closest_species_hit\tBLAST_description\n" %(datetime.date.today()))

    #call function to get the scaffold to gene dict
    scaffold_to_gene_dict, gene_to_exon_count, gene_start_stop_dict = parse_gff(gff)    
    #call function to get gene_set, gene_to_comment_dict
    HGT_predicted_gene_set, gene_to_comment_dict,\
        gene_to_HGTspeces_discription_dict, gene_to_AI, \
        gene_to_HGT_percent_identity = LTG_file(LTG)
    #get_scaffold_coverage from import coverage.py
    if bam_file:
        overall_coverage = "overall_coverage.txt"
        overall_expression_dic = get_total_coverage(bam_file, overall_coverage)
        HGT_gene_to_genic_cov_dic, scaffold_mean_SD_cov_dict, \
        mean_genomic_cov, standard_dev_genomic_cov = get_scaffold_coverage(genome, \
                                                    scaffold_to_gene_dict, gene_start_stop_dict,\
                                                    bam_file, overall_expression_dic,\
                                                    HGT_predicted_gene_set)

    #print gene_to_exon_count
    #call function to get rna seq mapping TPM
    gene_to_expression = parse_rnaseq(rnaseq)

    #call function with DNA file
    gene_AT_cont_dic, the_mean, standard_dev = get_stats_on_AT_content(dna_file)
    print "the AVR AT = %f with SD %f " %(the_mean, standard_dev)
    
    for gene, comment in gene_to_comment_dict.items():
        if bam_file:
            genomic_cov_from_mean = how_many_sd_from_mean(mean_genomic_cov, \
                                                      standard_dev_genomic_cov, \
                                                 HGT_gene_to_genic_cov_dic[gene])
        else:
            mean_genomic_cov = 0
            standard_dev_genomic_cov=0
            genomic_cov_from_mean=0
        try:
            AI = gene_to_AI[gene]
        except:
            ValueError
            AI = "NA"
            
        HGT_info_dataformatted = check_HGT_AT_vs_global_AT(gene_AT_cont_dic, AI, the_mean,
                                                           standard_dev, gene, comment,
                                                           sd_numbers, gene_to_expression, \
                                                           gene_to_exon_count,
                                                           gene_to_HGTspeces_discription_dict, \
                                                           genomic_cov_from_mean,
                                                           gene_to_HGT_percent_identity)
        f_out.write(HGT_info_dataformatted)

    for scaffold, genes in scaffold_to_gene_dict.items():
        descrption_to_add = ""
        genes_string = ""
        AI_values = ""
        bad_contig = True
        for gene in genes:
            #print gene
            if gene not in HGT_predicted_gene_set:
                bad_contig = False

        if bad_contig == True:
            #genes_string = ""
            genes_on_scaffold = scaffold_to_gene_dict[scaffold]
            #print "genes_on_scaffold", genes_on_scaffold
            for member in genes_on_scaffold:
                genes_string = genes_string+" "+member
                try:
                    AI = gene_to_AI[gene]
                    # HGT looswe threshold of 20 set here.. for further examination
                    if int(AI) < 15:
                        bad_contig = False
                        continue
                except:
                    ValueError
                    bad_contig = False
                    AI = "NA"
                    continue
                AI_values = AI_values+" "+AI

            if bad_contig == True:
                print ("Bad scaffold = %s" %(scaffold))
                descrpt = gene_to_HGTspeces_discription_dict[gene]
                descrption_to_add = descrption_to_add+" "+descrpt
                data_formatted = "%s\tBad_scaffold\t%s\t%s\t%s\n" %(scaffold,
                                                                    genes_string,
                                                                    AI_values,
                                                                    descrption_to_add)
                out.write(data_formatted)
    out.close()


#################################################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

Tool to refine the HGT predicted gene based on RNAseq cov, genomic cov, exon number, percentage
identity to best non-metazoan hit and AT content that differes from normal.

$ python ~/misc_python/Lateral_gene_transfer_prediction_tool/check_contaminants_on_contigs.py
        --gff ../augustus.gff3 -LTG LTG_LGT_candifates.out (default)

python ~/misc_python/Lateral_gene_transfer_prediction_tool/check_contaminants_on_contigs.py
    --bam sorted.bam --gff augustus.gff3 --LTG LTG_LGT_candifates_AI_30plus.out -s 0
    -r Rp.nt.fasta_quant.sf -g Rp.v1_alt.fasta --dna Rp.nt.fasta -o test


Requires:
samtools 1.2 or later for Bam file
Biopython
NUmpy

1) GENRATE bam file with genomic reads mapped to it:
How ever you want to do it, but sort and index your bam file
transrate --assembly genome.fasta --left genomic_reads.r1.fq.gz --right genomic_reads.r2.fq.gz --threads 12

BAM file is not need and can be run without it.  = much faster!!

2) GFF3
You may have to tidy and sort your GFF to a GFF3. Use GenomeTools
.. convert augustus.gft to gff3
.. gt gtf_to_gff3 -o test.gff3 -tidy augustus.gtf
or
.. gt gff3 -sort -tidy augustus.gff > formatted.gff3

3) LTG_LGT_candifates_AI_30plus.out:
This is the ouput from the Lateral_gene_transfer prediction tool. Precurser to this script.

4) RNAseq_coverage:
Agin, however you want to generate it. e.g. 
transrate --assembly gene.cds --left rnaseq_r1.fq.gz --right rnaseq_r2.fq.gz --threads 12

5) Genome seq -g

6) cds of genes:
If you dont have it can use:
gffread *gff -g genome.fasta -x nt.fa -y aa.fa

BAD SCAFFOLDS??

The script will check to see if a contig is only made up of LTG/HGT predicted genes.
If so, then this contig is suspect
and therefore should be considered as contimination.
Users are encouraged to used Blobplots of the genome assemblies before they get to this point.

"""

parser = OptionParser(usage=usage)

parser.add_option("--gff", "--GFF", dest="gff", default="test.gff",
                  help="gff file for predicted genes. ",
                  metavar="FILE")
parser.add_option("--LTG", dest="LTG", default="test_LTG_LGT_candifates.out",
                  help="LTG outfile. This is the output generated "
                  " from the Lateral_gene_transfer_prediction_tool "
                  " (LTG_results.out_Alien_index.out)",
                  metavar="FILE")
parser.add_option("--dna", dest="dna", default=None,
                  help="predicted cds nucleotide genes for AT content stats ",
                  metavar="FILE")
parser.add_option("-g", "--genome", dest="genome", default=None,
                  help="genome.fasta ",
                  metavar="FILE")
parser.add_option("-s", dest="sd_numbers", default=0,
                  help="the number of stadard deviations away from the mean"
                  " for identifying genes "
                  " that differ from normal AT content. default=0")
parser.add_option("-r", "--rnaseq", dest="rnaseq", default=None,
                  help="RNAseq expression profile for genes. "
                  " in format # Name	Length	TPM	NumReads "
                  " standard Sailfish output.")
parser.add_option("-b", "--bam", dest="bam_file", default=False,
                  help="bam file (sorted, indexed for samtools) "
                  " with genomic reads mapped to geneome "
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
# Run as script
if __name__ == '__main__':
    # call the main function
    if not os.path.isfile(gff):
        print("sorry, couldn't open the file: " % gff)
        print ("current working directory is :", os.getcwd() + "\n")
        print ("files are :", [f for f in os.listdir('.')])
        sys_exit("\n\nInput gff file not found: %s" % gff)


    if not os.path.isfile(LTG):
        sys_exit("Input LTG file not found: %s" % LTG)
    if not os.path.isfile(dna):
        sys_exit("Input dna file not found: %s" % dna)
        
        
    check_scaffolds_for_only_HGT_genes(genome, gff,
                                       LTG, dna,
                                       sd_numbers,
                                       rnaseq, bam_file,
                                       out_file)


