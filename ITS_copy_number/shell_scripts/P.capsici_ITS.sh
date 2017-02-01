#!/bin/bash
# put out files in current directory
#$ -cwd
#$ -l hostname="n13*"
# Deliver notifications to the following address
# Send notifications when the job begins (b), ends (e), or is aborted (a)
#dollar -M email.address@somewhere .... currently doesnt work
#dollar -m a 
#Abort on any error,
#set -euo pipefail
#(and thus quit the script right away)

echo Running on $HOSTNAME
echo Current PATH is $PATH
#source ~/.bash_profile

cd ~/scratch/tree_health/ITS_ratio
export TMP=~/scratch/${USER}_${JOB_ID}

##################################################################################################################################################################
# THESE VARIABLE NEED TO BE FILLED IN BY USER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

species=Phytophthora_capsici

genome_prefix=GCA_000149735.1_ASM14973v1

genome_fasta=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000149735.1_ASM14973v1/GCA_000149735.1_ASM14973v1_genomic.fna.gz

genome_GFF=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000149735.1_ASM14973v1/GCA_000149735.1_ASM14973v1_genomic.gff.gz
#http://www.ncbi.nlm.nih.gov/bioproject/?term=ADVJ00000000
#read_1_link=NO DATA AVAILBLE

#read_2_link=NO DATA AVAILBLE

trimmomatic_path=~/Downloads/Trimmomatic-0.32

SRA_prefix=SRR610900

repository_path=~/misc_python/THAPBI/ITS_region_genomic_coverage

num_threads=12


# NOTHING TO BE FILLED IN BY USER FROM HERE!!!!
##################################################################################################################################################################

export species
export genome_prefix
export genome_fasta
export genome_GFF
export read_1_link
export read_2_link
export trimmomatic_path
export SRA_prefix
export repository_path
export num_threads


./ITS_genomic_coverage_ratio_finding.sh
