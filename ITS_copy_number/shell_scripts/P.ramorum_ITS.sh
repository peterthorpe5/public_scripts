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

species=Phytophthora_ramorum

genome_prefix=Phytophthora_ramorum.ASM14973v1.31

genome_fasta=ftp://ftp.ensemblgenomes.org/pub/protists/release-31/fasta/phytophthora_ramorum/dna/Phytophthora_ramorum.ASM14973v1.31.dna.genome.fa.gz

genome_GFF=ftp://ftp.ensemblgenomes.org/pub/protists/release-31/gff3/phytophthora_ramorum/Phytophthora_ramorum.ASM14973v1.31.gff3.gz

read_1_link=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR610/SRR610900/SRR610900_1.fastq.gz

read_2_link=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR610/SRR610900/SRR610900_2.fastq.gz

#read_1_b_link=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR610/SRR610900/SRR610900_1.fastq.gz

#read_2_b_link=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR610/SRR610900/SRR610900_2.fastq.gz

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
