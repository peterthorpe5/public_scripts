#!/usr/bin/env bash
set -e
# Title: pipline to identify Phy species by clustering with ITS regions
# Author:  Peter Thorpe
# (C) The James Hutton Institute 2016


####################################################################################
# THESE VARIABLE NEED TO BE FILLED IN BY USER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Name_of_project=testing

# read name prefix not actually needed.
read_name_prefix="M01157:20:000000000-D07KA:1:"

left_read_file=DNAMIX_S95_L001_R1_001.fastq.gz

right_read_file=DNAMIX_S95_L001_R2_001.fastq.gz

trimmomatic_path=~/Downloads/Trimmomatic-0.32

repository_path=~/misc_python/THAPBI/Phyt_ITS_identifying_pipeline

num_threads=4

working_directory_path=$HOME/misc_python/THAPBI/Phyt_ITS_identifying_pipeline/data


# NOTHING TO BE FILLED IN BY USER FROM HERE!!!!
######################################################################################
#not needed but, may use it as a configuration style script later
export Name_of_project
export left_read_file
export right_read_file
export trimmomatic_path
export repository_path
export num_threads
export working_directory_path
export read_name_prefix


cd ${working_directory_path}

${repository_path}/Phyt_ITS_identify_pipline.sh

