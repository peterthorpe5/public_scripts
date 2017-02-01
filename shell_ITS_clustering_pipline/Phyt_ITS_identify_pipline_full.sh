#!/usr/bin/env bash
set -e
# Title: pipline to identify Phy species by clustering with ITS regions
# Author:  Peter Thorpe
# (C) The James Hutton Institute 2016


##################################################################################################################################################################
# THESE VARIABLE NEED TO BE FILLED IN BY USER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Name_of_project=testing

read_name_prefix="M01157:20:000000000-D07KA:1:"

left_read_file=DNAMIX_S95_L001_R1_001.fastq.gz

right_read_file=DNAMIX_S95_L001_R2_001.fastq.gz

trimmomatic_path=~/Downloads/Trimmomatic-0.32

repository_path=~/misc_python/THAPBI/Phyt_ITS_identifying_pipeline

num_threads=4

working_directory_path=$HOME/misc_python/THAPBI/Phyt_ITS_identifying_pipeline/data


# NOTHING TO BE FILLED IN BY USER FROM HERE!!!!
##################################################################################################################################################################
#not needed but, may use it as a configuration style script later
export trimmomatic_path
export repository_path
export num_threads


cd ${working_directory_path}


# Create directory for output
#mkdir fastq-join_joined
mkdir ${Name_of_project}_results
mkdir ${Name_of_project}_outfiles
mkdir ${Name_of_project}_quality_control_files

# QC "forward" and "reverse" reads
echo "running fastqc on raw reads"
cmd_fastqc="fastqc ${left_read_file}" 
echo ${cmd_fastqc}
eval ${cmd_fastqc}
cmd_fastqc2="fastqc ${right_read_file}" 
echo ${cmd_fastqc2}
eval ${cmd_fastqc2}
echo "cmd_fastqc done"

#quality trim the reads
echo "Trimming:"
echo "have you added to the adapter seq to the triming database file?"
echo "do we want to keep the adapters in right until the end to track the seq origin?"
cmd_trimming="java -jar ${trimmomatic_path}/trimmomatic-0.32.jar PE -threads ${num_threads} -phred33 ${left_read_file} ${right_read_file} ${Name_of_project}_R1.fq.gz unpaired_R1.fq.gz ${Name_of_project}_R2.fq.gz unpaired_R2.fq.gz ILLUMINACLIP:${repository_path}/database_files/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75" 
echo ${cmd_trimming}
eval ${cmd_trimming}
echo "Trimming done"

#remove unpaired reads
echo "removing unpaired files."
rm unpaired*

# Conduct QC on trimmed reads
echo "running fastqc on trimmed reads"
cmd_fastqc="fastqc ${Name_of_project}_R1.fq.gz" 
echo ${cmd_fastqc}
eval ${cmd_fastqc}
cmd_fastqc2="fastqc ${Name_of_project}_R2.fq.gz" 
echo ${cmd_fastqc2}
eval ${cmd_fastqc2}
echo "cmd_fastqc done"

# Join trimmed paired-end read files together
# use a program called PEAR
echo "puttin the reads together:"
cmd_pear="pear -f ${Name_of_project}_R1.fq.gz -r ${Name_of_project}_R2.fq.gz -o ${working_directory_path}/${Name_of_project}_outfiles/${Name_of_project}_PEAR" 
echo ${cmd_pear}
eval ${cmd_pear}
echo "cmd_pear done"

# Move the joined read files to the output subdirectory
#mv ${Name_of_project}_PEAR.assembled.fastq ./${Name_of_project}_outfiles

# Change directory to the output subdirectory
#cd ${Name_of_project}_outfiles

# Conduct QC on joined reads
fastqc ./${Name_of_project}_outfiles/${Name_of_project}_PEAR.assembled.fastq

# move all the fastqc files to a folder
wait

mv *.zip ./${Name_of_project}_quality_control_files
mv *.html ./${Name_of_project}_quality_control_files

# Convert read format from FASTQ to FASTA
# convert_format comes from seq_crumbs: https://github.com/JoseBlanca/seq_crumbs
echo "converting the fq to fa file"
cmd_convert="convert_format -t fastq -o ./${Name_of_project}_outfiles/fastqjoin.join.fasta -f fasta ./${Name_of_project}_outfiles/${Name_of_project}_PEAR.assembled.fastq" 
echo ${cmd_convert}
eval ${cmd_convert}
echo "cmd_convert done"



echo "running swarm: swarm [OPTIONS] [filename]"
#https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf
# add --append-abundance 1 - added this 
cmd_swarm="swarm -t ${num_threads} --append-abundance 1 -d 1 -o swarm_clustering_on_data_with_no_database ./${Name_of_project}_outfiles/fastqjoin.join.fasta" 
echo ${cmd_swarm}
eval ${cmd_swarm}
echo "cmd_swarm done"


# now to check some real data
# cat the assembled ITS and the PhyDB ITS sequence
echo "cat the assembled ITS and the PhyDB ITS sequences"
cmd_cat="cat ${repository_path}/database_files/Phy_ITSregions_all_20160601.coded_name.fasta ${repository_path}/data/${Name_of_project}_outfiles/fastqjoin.join.fasta > ${repository_path}/data/temp_reads_plus_database.fasta" 
echo ${cmd_cat}
eval ${cmd_cat}

# to avoid parsing gzip file. Uncompress these to a temp file. 
zcat ${working_directory_path}/${left_read_file} > ${working_directory_path}/temp_not_trimmedr1.fq
zcat ${working_directory_path}/${right_read_file} > ${working_directory_path}/temp_not_trimmedr2.fq

#run swarm with the different D values
values="1 2 3"
for v in ${values}
do
	# run the clustering on the REAL dataset
	echo "running swarm at ${v} thresholds"
	swarm_command="swarm -t ${num_threads} --append-abundance 1 -d ${v} -o ${Name_of_project}_reads_cluseterd_with_phy_ITS_Swarmd${v}.out ${repository_path}/data/temp_reads_plus_database.fasta" 
	echo ${swarm_command}
	eval ${swarm_command}
	wait
	
	# python summarise clussters
	echo "python to summarise clusters"
	cmd_summarise_clustering="python ${repository_path}/Python_ITS_scripts/parse_clusters.py -i ${Name_of_project}_reads_cluseterd_with_phy_ITS_Swarmd${v}.out -d ${repository_path}/database_files/name_to_species_database.txt -o ${working_directory_path}/${Name_of_project}_results/${Name_of_project}_reads_cluseterd_with_phy_ITS_Swarmd${v}.RESULTS" 
	echo ${cmd_summarise_clustering}
	eval ${cmd_summarise_clustering}
	wait
	
	# python graphs the  clussters
	echo "python to graphically represent clusters"
	cmd_summarise_graphs="python ${repository_path}/Python_ITS_scripts/draw_bar_chart_of_clusters.py -i ${working_directory_path}/${Name_of_project}_results/${Name_of_project}_reads_cluseterd_with_phy_ITS_Swarmd${v}.RESULTS -o ${working_directory_path}/${Name_of_project}_results/${Name_of_project}_reads_clust_w_phy_ITS_Swarmd${v}.graph" 
	echo ${cmd_summarise_graphs}
	eval ${cmd_summarise_graphs}
	wait
	
	# python find Phy species identified
	echo "python to find Phy species identified - I am using the original reads, not the trimmed ones to get the barcodes."
	cmd_what_phy_species="python ${repository_path}/Python_ITS_scripts/get_results_from_cluster.py --left ${working_directory_path}/temp_not_trimmedr1.fq --right ${working_directory_path}/temp_not_trimmedr2.fq -i ${working_directory_path}/${Name_of_project}_results/${Name_of_project}_reads_cluseterd_with_phy_ITS_Swarmd${v}.RESULTS -o ${working_directory_path}/${Name_of_project}_results/${Name_of_project}_Phytophthora_species_identified_swarm_${v}.txt" 
	echo ${cmd_what_phy_species}
	eval ${cmd_what_phy_species}
	wait
done

rm ${working_directory_path}/temp_not_trimmedr1.fq
rm ${working_directory_path}/temp_not_trimmedr2.fq








