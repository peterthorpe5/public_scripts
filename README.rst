Peter Thorpe's public bioinformatics scripts.
=============================================
for most programs here, help can be accessed by asking for help at the command line:

Type ``python script_name.py -h`` for how to use them.

This is an ever growing repository of tools which I have used/ using for the various projects I am involved in. 

With in here you will find:


Alternative_to_ITS1_finding
===========================
This was an early draft to try an identify single copy common regions within genomes which primers could be desinged for as an alternative metabarcoding region to ITS1.
This is not under any further development. 

genomic_upstream_regions
===========================
This gets the upstream regions of a given gene set to help identify promoter regions. Used in: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0985-1
The genome of the yellow potato cyst nematode, Globodera rostochiensis, reveals insights into the basis of parasitism and virulence      

blast_output
============
Some scripts in here to identify top blast hits, filter BLAST hits based on a given phylum (or whatever tax id is given).

identifying_pallindrome
=======================
A collection of scripts to identify and plot pallindromes in a given geneome sequence. Please read the Word file in the folder if you want to know more.

reformat_fasta_hints_names_for_Braker
=====================================
Script to rename the fasta names and hints names for BRAKER gene prediction. For me, this failed if I did not do this myself.

convert_file_format
===================
A collection of scripts to convert file formats from one type to another.
    
ITS_copy_number
===============
A pipeline and collection of scripts to estimate the copy number of ITS1 regions (or any other given gene of interest) based on genomic read coverage
    
shell_ITS_clustering_pipline
============================
A metabarcoding clustering pipeline wrote in shell. This is a draft for the upcoming metapy.py pipeline (https://github.com/widdowquinn/THAPBI-pycits/tree/master). 

Diamond_BLAST_add_taxonomic_info
================================
Tool to post taxonomically annotate a DIAMOND blast output. 
 
Lateral_gene_transfer_prediction_tool
=====================================
Tool to predict horizontal or lateral gene transfer.

split_up_fasta_file_into_N_files
================================
Tool to split up a large fasta file in N smaller fasta files

domain_searching
================
Pipeline to identity and align domains of interest.

NGS
===
Tools for working with Illumina data

transposon_analysis
===================
Tools and pipelines for transposon analysis in genomes.

Fix_five_prime
==============
Tool to refine the 5 prime start codon after Transdecoder has predicted the CDS from an RNAseq assembly
 
primer_designer
===============
Under development

gene_model_testing
==================
Pipeline to tests gene models and gain information as to how good they are. They is no one method to do this!
      
produce_random_seq
===================
Program to produce N number of random sequences with the average length and average GC of a given database. 

