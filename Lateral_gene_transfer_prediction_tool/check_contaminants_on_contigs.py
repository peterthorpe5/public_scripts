#!/usr/bin/env python
#author: Peter Thorpe September 2015. The James Hutton Insitute, Dundee, UK.

#imports
import os
import sys
from sys import stdin,argv
import sys
from optparse import OptionParser

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


def parse_gff(gff):
    "function to parse GFF and produce a scaffold_to_gene_dict"
    f_in = open(gff, "r")
    # data to output of function
    scaffold_to_gene_dict = dict()
    #iterate through gff
    for line in f_in:
        if line.startswith("#"):
            continue
        if not line.split("\t")[2] == "gene":
            continue
        scaffold,a,b,c,d,e,f,g,gene_info = line.split("\t")
        gene = gene_info.replace(";", "").split("ID=")[1]
        gene = gene.split(".gene")[0]
        gene = gene.rstrip("\n")
        #scaffold_to_gene_dict[scaffold]=gene.rstrip("\n")
        if not scaffold in scaffold_to_gene_dict:
            scaffold_to_gene_dict[scaffold]=[gene]
        else:
            scaffold_to_gene_dict[scaffold].append(gene)
    #print scaffold_to_gene_dict
    f_in.close()
    return scaffold_to_gene_dict

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
        gene = gene.replace("ID=", "").split("gene=g")[0]
        gene = gene.rstrip()
        if not ".t" in gene:
            gene = gene.replace("t", ".t")
        HGT_predicted_gene_set.add(gene)
        gene_to_comment_dict[gene] = comment.rstrip("\n")
    #print HGT_predicted_gene_set
    f_in.close()
    return HGT_predicted_gene_set, gene_to_comment_dict

# main function

def check_scaffolds_for_only_HGT_genes(gff, LTG, out_file):
    """main function. This calls the other function to get a dictionary
    of genes on scaffolds. A list of HGT/LTG genes and check the scaffolds
    to identify those that only have HGT genes on them. If so, then this
    is most likely a contaminant contig/scaffold"""
    out = open(out_file, "w")
    
    #call function to get the scaffold to gene dict
    scaffold_to_gene_dict = parse_gff(gff)
    #call function to get gene_set, gene_to_comment_dict
    HGT_predicted_gene_set, gene_to_comment_dict = LTG_file(LTG)
    #for scaffold in scaffold_to_gene_dict:
        #bad_contig = True
        #for gene in scaffold:
            #if gene not in HGT_predicted_gene_set:
                #bad_contig = False
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



You may have to tidy and sort your GFF to a GFF3. Use GenomeTools

STEPS 1)

convert augustus.gft to gff3

gt gtf_to_gff3 -o test.gff3 -tidy augustus.gtf


or


gt gff3 -sort -tidy augustus.gff > formatted.gff3

. 

"""

parser = OptionParser(usage=usage)

parser.add_option("--gff", dest="gff", default="test.gff",
                  help="hintsfile",
                  metavar="FILE")
parser.add_option("--LTG", dest="LTG", default="test_LTG_LGT_candifates.out",
                  help="LTG outfile. ",
                  metavar="FILE")

parser.add_option("-o", "--out_file", dest="out_file", default="test_output.txt",
                  help="outfile to list the bad contigs")


(options, args) = parser.parse_args()


gff = options.gff
out_file = options.out_file
LTG = options.LTG


#run the program

if not os.path.isfile(gff):
    sys_exit("Input gff file not found: %s" % gff)

if not os.path.isfile(LTG):
    sys_exit("Input LTG file not found: %s" % LTG)
    
check_scaffolds_for_only_HGT_genes(gff, LTG, out_file)


