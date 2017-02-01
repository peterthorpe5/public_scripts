#!/usr/bin/env python
#title: lazy mans way of generating a load of config files
#in the cluster
#author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

usage  = """ THIS is a python 2.7 script which will make a load of .sh scripts given
a
 list of read file. 

USER must copy and paste a list of read_file_names into this file
after the 'list_of_read_files = quotequotequote variable.

If your reads do not end in 001, please change lines:

 >prefix = i.split("001.f")[0]
and
> temp_name = i.split("_L001_R1_")[0]+".sh"

"""
print usage

# user paste read_files here ... keep three quote before and after please.
list_of_read_files = """
DNAMIX_S95_L001_R1_001.fastq.gz      GSB-301111_S65_L001_R1_001.fastq.gz  IGB-180712_S17_L001_R1_001.fastq.gz
DNAMIX_S95_L001_R2_001.fastq.gz      GSB-301111_S65_L001_R2_001.fastq.gz  IGB-180712_S17_L001_R2_001.fastq.gz
ECN_S1_L001_R1_001.fastq.gz          IGB-010212_S5_L001_R1_001.fastq.gz   IGB-180713_S37_L001_R1_001.fastq.gz
ECN_S1_L001_R2_001.fastq.gz          IGB-010212_S5_L001_R2_001.fastq.gz   IGB-180713_S37_L001_R2_001.fastq.gz
GSB-010513_S80_L001_R1_001.fastq.gz  IGB-010413_S29_L001_R1_001.fastq.gz  IGB-190813_S39_L001_R1_001.fastq.gz
GSB-010513_S80_L001_R2_001.fastq.gz  IGB-010413_S29_L001_R2_001.fastq.gz  IGB-190813_S39_L001_R2_001.fastq.gz
GSB-020512_S73_L001_R1_001.fastq.gz  IGB-010813_S38_L001_R1_001.fastq.gz  IGB-200112_S4_L001_R1_001.fastq.gz
GSB-020512_S73_L001_R2_001.fastq.gz  IGB-010813_S38_L001_R2_001.fastq.gz  IGB-200112_S4_L001_R2_001.fastq.gz
GSB-021013_S87_L001_R1_001.fastq.gz  IGB-030113_S24_L001_R1_001.fastq.gz  IGB-200612_S15_L001_R1_001.fastq.gz
GSB-021013_S87_L001_R2_001.fastq.gz  IGB-030113_S24_L001_R2_001.fastq.gz  IGB-200612_S15_L001_R2_001.fastq.gz
GSB-040412_S71_L001_R1_001.fastq.gz  IGB-030214_S49_L001_R1_001.fastq.gz  IGB-220512_S13_L001_R1_001.fastq.gz
GSB-040412_S71_L001_R2_001.fastq.gz  IGB-030214_S49_L001_R2_001.fastq.gz  IGB-220512_S13_L001_R2_001.fastq.gz
GSB-040913_S86_L001_R1_001.fastq.gz  IGB-030812_S18_L001_R1_001.fastq.gz  IGB-231112_S21_L001_R1_001.fastq.gz
GSB-040913_S86_L001_R2_001.fastq.gz  IGB-030812_S18_L001_R2_001.fastq.gz  IGB-231112_S21_L001_R2_001.fastq.gz
GSB-050214_S92_L001_R1_001.fastq.gz  IGB-031012_S20_L001_R1_001.fastq.gz  IGB-250613_S35_L001_R1_001.fastq.gz
GSB-050214_S92_L001_R2_001.fastq.gz  IGB-031012_S20_L001_R2_001.fastq.gz  IGB-250613_S35_L001_R2_001.fastq.gz
GSB-050314_S94_L001_R1_001.fastq.gz  IGB-040213_S25_L001_R1_001.fastq.gz  IGB-251113_S45_L001_R1_001.fastq.gz
GSB-050314_S94_L001_R2_001.fastq.gz  IGB-040213_S25_L001_R2_001.fastq.gz  IGB-251113_S45_L001_R2_001.fastq.gz
GSB-060313_S78_L001_R1_001.fastq.gz  IGB-040313_S27_L001_R1_001.fastq.gz  IGB-260412_S11_L001_R1_001.fastq.gz
GSB-060313_S78_L001_R2_001.fastq.gz  IGB-040313_S27_L001_R2_001.fastq.gz  IGB-260412_S11_L001_R2_001.fastq.gz
GSB-070312_S69_L001_R1_001.fastq.gz  IGB-050712_S16_L001_R1_001.fastq.gz  IGB-270513_S33_L001_R1_001.fastq.gz
GSB-070312_S69_L001_R2_001.fastq.gz  IGB-050712_S16_L001_R2_001.fastq.gz  IGB-270513_S33_L001_R2_001.fastq.gz
GSB-100713_S83_L001_R1_001.fastq.gz  IGB-060911_S2_L001_R1_001.fastq.gz   IGB-280312_S9_L001_R1_001.fastq.gz
GSB-100713_S83_L001_R2_001.fastq.gz  IGB-060911_S2_L001_R2_001.fastq.gz   IGB-280312_S9_L001_R2_001.fastq.gz
GSB-110112_S67_L001_R1_001.fastq.gz  IGB-070512_S12_L001_R1_001.fastq.gz  IGB-281013_S43_L001_R1_001.fastq.gz
GSB-110112_S67_L001_R2_001.fastq.gz  IGB-070512_S12_L001_R2_001.fastq.gz  IGB-281013_S43_L001_R2_001.fastq.gz
GSB-111213_S91_L001_R1_001.fastq.gz  IGB-071212_S22_L001_R1_001.fastq.gz  IGB-290212_S7_L001_R1_001.fastq.gz
GSB-111213_S91_L001_R2_001.fastq.gz  IGB-071212_S22_L001_R2_001.fastq.gz  IGB-290212_S7_L001_R2_001.fastq.gz
GSB-120613_S82_L001_R1_001.fastq.gz  IGB-080114_S48_L001_R1_001.fastq.gz  IGB-290413_S31_L001_R1_001.fastq.gz
GSB-120613_S82_L001_R2_001.fastq.gz  IGB-080114_S48_L001_R2_001.fastq.gz  IGB-290413_S31_L001_R2_001.fastq.gz
GSB-130612_S75_L001_R1_001.fastq.gz  IGB-080612_S14_L001_R1_001.fastq.gz  IGB-300813_S40_L001_R1_001.fastq.gz
GSB-130612_S75_L001_R2_001.fastq.gz  IGB-080612_S14_L001_R2_001.fastq.gz  IGB-300813_S40_L001_R2_001.fastq.gz
GSB-141211_S66_L001_R1_001.fastq.gz  IGB-080713_S36_L001_R1_001.fastq.gz  IGB-300913_S42_L001_R1_001.fastq.gz
GSB-141211_S66_L001_R2_001.fastq.gz  IGB-080713_S36_L001_R2_001.fastq.gz  IGB-300913_S42_L001_R2_001.fastq.gz
GSB-150513_S81_L001_R1_001.fastq.gz  IGB-091213_S46_L001_R1_001.fastq.gz  PCRMIX_S96_L001_R1_001.fastq.gz
GSB-150513_S81_L001_R2_001.fastq.gz  IGB-091213_S46_L001_R2_001.fastq.gz  PCRMIX_S96_L001_R2_001.fastq.gz
GSB-160512_S74_L001_R1_001.fastq.gz  IGB-100613_S34_L001_R1_001.fastq.gz  SRB-030713_S56_L001_R1_001.fastq.gz
GSB-160512_S74_L001_R2_001.fastq.gz  IGB-100613_S34_L001_R2_001.fastq.gz  SRB-030713_S56_L001_R2_001.fastq.gz
GSB-161013_S88_L001_R1_001.fastq.gz  IGB-110112_S3_L001_R1_001.fastq.gz   SRB-050613_S55_L001_R1_001.fastq.gz
GSB-161013_S88_L001_R2_001.fastq.gz  IGB-110112_S3_L001_R2_001.fastq.gz   SRB-050613_S55_L001_R2_001.fastq.gz
GSB-161111_S64_L001_R1_001.fastq.gz  IGB-111113_S44_L001_R1_001.fastq.gz  SRB-120214_S62_L001_R1_001.fastq.gz
GSB-161111_S64_L001_R2_001.fastq.gz  IGB-111113_S44_L001_R2_001.fastq.gz  SRB-120214_S62_L001_R2_001.fastq.gz
GSB-170413_S79_L001_R1_001.fastq.gz  IGB-120412_S10_L001_R1_001.fastq.gz  SRB-120314_S63_L001_R1_001.fastq.gz
GSB-170413_S79_L001_R2_001.fastq.gz  IGB-120412_S10_L001_R2_001.fastq.gz  SRB-120314_S63_L001_R2_001.fastq.gz
GSB-180412_S72_L001_R1_001.fastq.gz  IGB-130513_S32_L001_R1_001.fastq.gz  SRB-130313_S53_L001_R1_001.fastq.gz
GSB-180412_S72_L001_R2_001.fastq.gz  IGB-130513_S32_L001_R2_001.fastq.gz  SRB-130313_S53_L001_R2_001.fastq.gz
GSB-190214_S93_L001_R1_001.fastq.gz  IGB-140312_S8_L001_R1_001.fastq.gz   SRB-150114_S61_L001_R1_001.fastq.gz
GSB-190214_S93_L001_R2_001.fastq.gz  IGB-140312_S8_L001_R2_001.fastq.gz   SRB-150114_S61_L001_R2_001.fastq.gz
GSB-210312_S70_L001_R1_001.fastq.gz  IGB-150212_S6_L001_R1_001.fastq.gz   SRB-180712_S52_L001_R1_001.fastq.gz
GSB-210312_S70_L001_R2_001.fastq.gz  IGB-150212_S6_L001_R2_001.fastq.gz   SRB-180712_S52_L001_R2_001.fastq.gz
GSB-210813_S85_L001_R1_001.fastq.gz  IGB-150413_S30_L001_R1_001.fastq.gz  SRB-230512_S51_L001_R1_001.fastq.gz
GSB-210813_S85_L001_R2_001.fastq.gz  IGB-150413_S30_L001_R2_001.fastq.gz  SRB-230512_S51_L001_R2_001.fastq.gz
GSB-240713_S84_L001_R1_001.fastq.gz  IGB-150812_S19_L001_R1_001.fastq.gz  SRB-231013_S60_L001_R1_001.fastq.gz
GSB-240713_S84_L001_R2_001.fastq.gz  IGB-150812_S19_L001_R2_001.fastq.gz  SRB-231013_S60_L001_R2_001.fastq.gz
GSB-250112_S68_L001_R1_001.fastq.gz  IGB-160913_S41_L001_R1_001.fastq.gz  SRB-250913_S59_L001_R1_001.fastq.gz
GSB-250112_S68_L001_R2_001.fastq.gz  IGB-160913_S41_L001_R2_001.fastq.gz  SRB-250913_S59_L001_R2_001.fastq.gz
GSB-250712_S76_L001_R1_001.fastq.gz  IGB-161213_S47_L001_R1_001.fastq.gz  SRB-280312_S50_L001_R1_001.fastq.gz
GSB-250712_S76_L001_R2_001.fastq.gz  IGB-161213_S47_L001_R2_001.fastq.gz  SRB-280312_S50_L001_R2_001.fastq.gz
GSB-271113_S90_L001_R1_001.fastq.gz  IGB-171212_S23_L001_R1_001.fastq.gz  SRB-280813_S58_L001_R1_001.fastq.gz
GSB-271113_S90_L001_R2_001.fastq.gz  IGB-171212_S23_L001_R2_001.fastq.gz  SRB-280813_S58_L001_R2_001.fastq.gz
GSB-281112_S77_L001_R1_001.fastq.gz  IGB-180213_S26_L001_R1_001.fastq.gz  SRB-290513_S54_L001_R1_001.fastq.gz
GSB-281112_S77_L001_R2_001.fastq.gz  IGB-180213_S26_L001_R2_001.fastq.gz  SRB-290513_S54_L001_R2_001.fastq.gz
GSB-301013_S89_L001_R1_001.fastq.gz  IGB-180313_S28_L001_R1_001.fastq.gz  SRB-310713_S57_L001_R1_001.fastq.gz
GSB-301013_S89_L001_R2_001.fastq.gz  IGB-180313_S28_L001_R2_001.fastq.gz  SRB-310713_S57_L001_R2_001.fastq.gz""".split()


seen_set = set([])


for i in sorted(list_of_read_files):
	prefix = i.split("001.f")[0]
	if prefix not in seen_set:
		#print prefix
		seen_set.add(prefix)

for i in sorted(seen_set):
	if "R1" in i:
		temp_name = i.split("_L001_R1_")[0]+".sh"
		f =open(temp_name, "w")
		print >> f, """#!/usr/bin/env bash
set -e
# Title: pipline to identify Phy species by clustering with ITS regions
# Author:  Peter Thorpe
# (C) The James Hutton Institute 2016


####################################################################################
# THESE VARIABLE NEED TO BE FILLED IN BY USER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"""
		print >> f, "Name_of_project='%s'" %(temp_name.replace("-", "_").split(".sh")[0])
		print >> f, "Name_of_ITS_database='ITS_database_NOT_confirmed_correct_last14bases_removed.fasta'"
		print >> f, "#Name_of_ITS_database='Phy_ITSregions_all_20160601.coded_name.fasta'"

		print >> f, "\nleft_read_file=$HOME/scratch/tree_health/201606_THAPBI_MiSeq_result/%s001.fastq.gz" %(i)
	if "R2" in i:
		print >> f,"\nright_read_file=$HOME/scratch/tree_health/201606_THAPBI_MiSeq_result/%s001.fastq.gz\n" %(i)
		print >> f, """trimmomatic_path=~/Downloads/Trimmomatic-0.32

repository_path=~/misc_python/THAPBI/Phyt_ITS_identifying_pipeline

num_threads=4

values="1 2 3"

# 53, not sure if this is correct. 
left_primer_length=53
right_primer_length=0
# need to check to see if the barcode len needs to be added as well as the primer len
barcode_length=0
read_prefix=M01157

working_directory_path=$HOME/scratch/tree_health/testing_pipline/%s

path_to_spades=$HOME/scratch/Downloads/SPAdes-3.9.0-Linux/bin/

# NOTHING TO BE FILLED IN BY USER FROM HERE!!!!
######################################################################################
#not needed but, may use it as a configuration style script later
export Name_of_project
export Name_of_ITS_database
export left_read_file
export right_read_file
export trimmomatic_path
export repository_path
export num_threads
export working_directory_path
export values
export left_primer_length
export right_primer_length 
export barcode_length
export read_prefix
export path_to_spades
#export read_name_prefix

mkdir ${working_directory_path}

cd ${working_directory_path}

# If you want error correction uncomment this one and comment out the following one
#${repository_path}/Phyt_ITS_identify_pipline_with_error_correction.sh

${repository_path}/Phyt_ITS_identify_pipline.sh


""" % (temp_name.replace("-", "_").split(".sh")[0])
		f.close()
