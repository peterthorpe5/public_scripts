READme: Sanger_read_metagenetics
================================

basic usage:

Go into a directory which has a load of .abi or .ab1 files which need to be classifed against a Phytophthora database.
Add metapy_automatic.py to your path, along with metapy_sanger_read.py (If you do this you will have to modify
the metapy_automatic.py code to suit your system and path). Simply run metapy_automatic.py

``./metapy_sanger_read.py`` -h 

metapy_automatic.py: Pycits/ metapy classify OTU using Sanger ab1 files: v1.0.0.
This script will run the method of all .abi or .ab1 files in a folder.
This come with a Phytophora database, which are all the enteries in 
NCBI as of Jan 2017. If you want your own database you will have to run
metapy_sanger_read.py -d your_datase.fasta


author: Peter Thorpe and Leighton Pritchard September 2017 . The James Hutton Insitute, Dundee, UK.

usage: metapy_sanger_read.py [--thread THREADS] [-a AB1] [-d OTU_DB]
                             [-e EVALUE] [-m MISMATCHES]
                             [--left_trim LEFT_TRIM] [--right_trim RIGHT_TRIM]
                             [--adaptors ADAPTORS] [--phred PHRED] [--verbose]
                             [--align] [--muscle MUSCLE]
                             [--trimmomatic TRIMMOMATIC] [--logfile LOGFILE]
                             [--cleanup] [-h] [--version]

Pipeline: cluster data for metabarcoding

optional arguments:
  --thread THREADS      number of threads
  -a AB1, --ab1 AB1     sanger ab1 file
  -d OTU_DB, --OTU_DB OTU_DB
                        database of seq of to compare against
  -e EVALUE, --evalue EVALUE
                        evalue to filter results with
  -m MISMATCHES, --mismatches MISMATCHES
                        number of mismatches to filter results with
  --left_trim LEFT_TRIM
                        left_trim for primers or conserved regions. Default 25
  --right_trim RIGHT_TRIM
                        right_trim for primers or conserved regions. Default 0
  --adaptors ADAPTORS   adaptors for trimming. Can supply custom file if
                        desired
  --phred PHRED         phred33 is default. Dont change unless sure
  --verbose             Report verbose output
  --align               to align clusters in the output you must have muscle
                        in your PATH as muscle
  --muscle MUSCLE       Path to MUSCLE... If version alreadyin PATH then leave
                        blank
  --trimmomatic TRIMMOMATIC
                        Path to trimmomatic... If version alreadyin PATH then
                        leave blank
  --logfile LOGFILE     Logfile name
  --cleanup             deletes most files the program creates
  -h, --help            Displays this help message type --version for version
  --version             show program's version number and exit
