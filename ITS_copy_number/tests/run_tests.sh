#!/bin/bash
# put out files in current directory
#$ -cwd
#$ -l hostname="n13*"
# Deliver notifications to the following address
# Send notifications when the job begins (b), ends (e), or is aborted (a)
#Abort on any error,
set -e
#(and thus quit the script right away)

#source ~/.bash_profile
#export PATH=$HOME/bin:$PATH
#echo Revised PATH is $PATH


cd ${HOME}/misc_python/THAPBI/ITS_region_genomic_coverage/tests

python ../filter_GFF.py --gff test_get_representative_blast_GFF_region.gff -o CURRENT_TEST.gff


if diff CURRENT_TEST.gff test_get_representative_blast_GFF_region_RESULTS.gff ; then echo "Differences found failed"; false; fi