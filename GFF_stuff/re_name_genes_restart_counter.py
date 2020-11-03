#!/usr/bin/env python
# Title:
# script to reformt gene names after gene prediction.
# GT sometime returns . as a name for all genes.
import sys
from optparse import OptionParser
import re


def reformat_gff_column(gff):
    """function to rename the gene in a GFF file. The automatically
    gnerated Augustus names/Braker name need to be altered.
    These are altered to perfix00001 - then 5 interger number gene name"""
    f_in = open(gff, "r")
    outfile = "fixed.gff"
    f_out = open(outfile, 'w')
    gene_count = 0
    prefix =  ""
    for line in f_in:
        if line.startswith("#"):
            f_out.write(line.rstrip() + "\n")
            continue
        #line = line.split("product=")[0]
        #line = line.split("Name=")[0]
        # split the coloumn up in the gff based on \t
        scaf, aug, info, start, stop, stats, direction,\
                    more_info, gene_info = line.split("\t")
        gene_info = gene_info.rstrip()
        
        if info == "gene":
            old_gene_num = gene_info.split(";Name=")[0]
            old_gene_num = old_gene_num.rstrip()
            old_gene_num = old_gene_num.split("GPALLIND_")[1]
            old_gene_num = old_gene_num.replace(";", "")
            
        #count the genes
        if info == "gene":
            
            gene_count = gene_count + 1
            gene_number = "%s%05d" %(prefix, gene_count)
        #gene_number = "%s%05d" %(prefix, gene_count)
        # old_names = re.compile("g[0-9]+\")
        #print("gene_info_line = ", gene_info)
        #print("genenum = ", old_gene_num, " replacing with ", gene_number)
        gene_info = gene_info.replace(old_gene_num, gene_number)
         
        data = "\t".join([scaf,
                          aug,
                          info,
                          start,
                          stop,stats,\
                          direction,
                          more_info,\
                          gene_info + "\n"])
        f_out.write(data)
    f_out.close()
    f_in.close()

reformat_gff_column("Globodera_pallida_LINDLEY.gff3")
