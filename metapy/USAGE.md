Usage instruction for metapy.py
===============================

It is a python3 tool and because it runs lots of tools it has a lot of requirements. You dont need to install them all.
Minimum you need ``trimmomatic``, ``Swarm``, ``pear``, ``biopython``

Run this to see the command options:
====================================
``metapy.py`` -h

usage: metapy.py 

Pipeline: cluster data for metabarcoding

options:
  --thread THREADS      number of threads
  -l LEFT, --left LEFT  left illumina reads, default is for tests
  -r RIGHT, --right RIGHT
                        right illumina reads
  -d OTU_DB, --OTU_DB OTU_DB
                        right illumina reads
  -a {pear,flash}, --assemble {pear,flash}
                        program to assemble with flash or pear default: pear
  --adaptors ADAPTORS   adaptors for trimming. Can supply custom file if
                        desired
  --left_trim LEFT_TRIM
                        left_trim for primers or conserved regions. Default 53
  --right_trim RIGHT_TRIM
                        right_trim for primers or conserved regions. Default
                        20 (5.8 region)
  --phred PHRED         phred33 is default. Dont change unless sure
  --cdhit_threshold CDHIT_THRESHOLD
                        percentage identify for cd-hit Default -0.99
  --swarm_d_value SWARM_D_VALUE
                        the difference d value for clustering in swarm.
                        Default 1
  --blastclust_threshold BLASTCLUST_THRESHOLD
                        the threshold for blastclust clustering Default -S
                        0.90
  --vesearch_threshold VESEARCH_THRESHOLD
                        the threshold for vsearch clustering Default 0.99
  --verbose             Report verbose output
  -e, --error_correction
                        to perform Illumina error correction
  --align               to align clusters in the output you must have muscle
                        in your PATH as muscle
  --percent_identity    blast the cluster to return pairwise percentage
                        identity
  --min_novel_cluster_threshold MIN_NOVEL_CLUSTER_THRESHOLD
                        min size of a cluster to consider as real anything
                        smaller than this is ignored
  --blastclust BLASTCLUST
                        Path to blastclust ... If version alreadyin PATH then
                        leave blank
  --muscle MUSCLE       Path to MUSCLE... If version alreadyin PATH then leave
                        blank
  --flash FLASH         Path to flash... If version alreadyin PATH then leave
                        blank
  --pear PEAR           Path to pear... If version alreadyin PATH then leave
                        blank
  --cd-hit-est CD_HIT   Path to cd-hit-est... If version alreadyin PATH then
                        leave blank
  --bowtie2 BOWTIE2     Path to bowtie2... If version alreadyin PATH then
                        leave blank
  --fastqc FASTQC       Path to fastqc... If version alreadyin PATH then leave
                        blank
  --spades SPADES       Path to spades.py... If version alreadyin PATH then
                        leave blank
  --vsearch VSEARCH     Path to vsearch... If version alreadyin PATH then
                        leave blank
  --trimmomatic TRIMMOMATIC
                        Path to trimmomatic... If version alreadyin PATH then
                        leave blank
  --swarm SWARM         Path to swarm... If version alreadyin PATH then leave
                        blank
  --samtools SAMTOOLS   Path to samtools... If version alreadyin PATH then
                        leave blank
  --logfile LOGFILE     Logfile name
  --Run_blastclust RUN_BLASTCLUST
                        Run_blastclust yes or no, This is slow default no
  --Run_cdhit RUN_CDHIT
                        Run_cdhit yes or no, This is slow default no
  --Run_dada2 RUN_DADA2
                        Run_dada2 yes or no, default no
  --plot PLOTS          create graphs of clusters yes or no, This is slow.
                        Default no.
  --pvalue PVALUE       pvalue for comparing the database lengths versus the
                        assembled lengths. At which point are the different?
  --standard_deviation STD
                        standard_deviation threshold for comparing the
                        assembled size versus the database sequence sizes.
                        This is to check database is sensible for the data
                        input
  --cleanup             deletes most files the program creates
  --qc                  performs QC at various stages. Turn off by: --qc False
  -h, --help            Displays this help message type --version for version
  --version             show program's version number and exit
