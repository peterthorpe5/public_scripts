#!/usr/bin/env bash
set -e
# Title: pipline to identify Phy species by clustering with ITS regions
# Author:  Peter Thorpe
# (C) The James Hutton Institute 2016


##################################################################################################################################################################
# THESE VARIABLE NEED TO BE FILLED IN BY USER  in the config file


#info for user:
# this is the main script. Nothing should be changed in here. 
# These variables are called by an individaul script see for and
#example Phy_CONFIG_file.sh


cd ${working_directory_path}

echo Name_of_project = ${Name_of_project}
echo Name_of_ITS_database = ${Name_of_ITS_database}
echo left_read_file = ${left_read_file}
echo right_read_file = ${right_read_file}
echo trimmomatic_path = ${trimmomatic_path}
echo repository_path = ${repository_path}
echo num_threads = ${num_threads}
echo working_directory_path = ${working_directory_path}
echo values = ${values}
echo left_primer_length = ${left_primer_length}
echo right_primer_length  = ${right_primer_length}
echo barcode_length = ${barcode_length}
echo read_prefix = ${read_prefix}
echo path_to_spades= ${path_to_spades}


# Create directory for output
#mkdir fastq-join_joined
echo "making dir ${Name_of_project}_results"
mkdir ${Name_of_project}_results
mkdir ${Name_of_project}_outfiles
mkdir ${Name_of_project}_quality_control_files

# QC "forward" and "reverse" reads
echo "1) running fastqc on raw reads"
cmd_fastqc="fastqc ${left_read_file}" 
echo ${cmd_fastqc}
eval ${cmd_fastqc}
cmd_fastqc2="fastqc ${right_read_file}" 
echo ${cmd_fastqc2}
eval ${cmd_fastqc2}
echo "cmd_fastqc done"

#quality trim the reads
echo "2) Trimming:"
echo "have you added to the adapter seq to the triming database file?"
echo "do we want to keep the adapters in right until the end to track the seq origin?"
cmd_trimming="java -jar ${trimmomatic_path}/trimmomatic-0.32.jar PE -threads ${num_threads} 
-phred33 ${left_read_file} ${right_read_file} ${Name_of_project}_R1.fq.gz unpaired_R1.fq.gz ${Name_of_project}_R2.fq.gz 
unpaired_R2.fq.gz 
ILLUMINACLIP:${repository_path}/database_files/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75" 
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

#error correction on reads using SPAdes
echo "2b) Bayes hammer error correction. Module in spades. "
error_correct="${path_to_spades}/spades.py --only-error-correction 
-1 ${Name_of_project}_R1.fq.gz 
-2 ${Name_of_project}_R2.fq.gz 
-o ${working_directory_path}/${Name_of_project}_error_correction"
echo ${error_correct}
eval ${error_correct}

mv ${working_directory_path}/${Name_of_project}_error_correction/corrected/${Name_of_project}* ${working_directory_path}
rm -rf ${working_directory_path}/${Name_of_project}_error_correction
rm *unpaired*

# Join trimmed paired-end read files together
# use a program called PEAR
echo "3) puttin the error corrected reads together:"
cmd_pear="pear -f ${Name_of_project}_R1.fq.00.0_0.cor.fastq.gz 
-r ${Name_of_project}_R2.fq.00.0_0.cor.fastq.gz 
-o ${working_directory_path}/${Name_of_project}_outfiles/${Name_of_project}_PEAR" 
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
echo "4) converting the fq to fa file"
cmd_convert="convert_format -t fastq 
-o ./${Name_of_project}_outfiles/fastqjoin.join.fasta_with_barcodes 
-f fasta ./${Name_of_project}_outfiles/${Name_of_project}_PEAR.assembled.fastq" 
echo ${cmd_convert}
eval ${cmd_convert}
echo "cmd_convert done"

# python remove barcodes and excess seq from sequences
echo "5) python to remove barcodes and excess seq"
cmd_remove_bar_excess="python ${repository_path}/Python_ITS_scripts/remove_barcodes.py 
-f ./${Name_of_project}_outfiles/fastqjoin.join.fasta_with_barcodes -l ${left_primer_length}
 -r ${right_primer_length} 
 --barcode ${barcode_length} 
 --ITS CCACAC 
 -o ./${Name_of_project}_outfiles/temp" 
echo ${cmd_remove_bar_excess}
eval ${cmd_remove_bar_excess}
wait

# python remove barcodes and excess seq from sequences
echo "5b) python rewite as fasta"
cmd_rewite="python ${repository_path}/Python_ITS_scripts/rewrite_as_fastafile.py 
-i ./${Name_of_project}_outfiles/temp 
-o ./${Name_of_project}_outfiles/fastqjoin.join.fasta" 
echo ${cmd_rewite}
eval ${cmd_rewite}
wait

echo "6) running swarm: swarm [OPTIONS] [filename]"
#https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf
# add --append-abundance 1 - added this 
cmd_swarm="swarm -t ${num_threads} --append-abundance 1 -d 1 
-o swarm_clustering_on_data_with_no_database 
./${Name_of_project}_outfiles/fastqjoin.join.fasta" 
echo ${cmd_swarm}
eval ${cmd_swarm}
echo "cmd_swarm done"


# now to check some real data
# cat the assembled ITS and the PhyDB ITS sequence
echo "7) cat the assembled ITS and the PhyDB ITS sequences"
cmd_cat="cat ${repository_path}/database_files/${Name_of_ITS_database} 
${working_directory_path}/${Name_of_project}_outfiles/fastqjoin.join.fasta 
> ${working_directory_path}/temp_reads_plus_database_full_names.fasta" 
echo ${cmd_cat}
eval ${cmd_cat}

# python rename the files so the names work in Swarm
echo "8) python rewite as fasta with coded names"
cmd_rename="python ${repository_path}/Python_ITS_scripts/completely_rename_ID.py 
-f ${working_directory_path}/temp_reads_plus_database_full_names.fasta 
-d ./${Name_of_project}_outfiles/databse_old_names_to_temp.txt 
-o ${working_directory_path}/temp_reads_plus_database.fasta" 
echo ${cmd_rename}
eval ${cmd_rename}
wait

# cat the assembled ITS and the PhyDB ITS database
echo "8b) cat the assembled ITS and the PhyDB ITS databases"
cmd_cat="cat ./${Name_of_project}_outfiles/databse_old_names_to_temp.txt 
${repository_path}/database_files/name_to_species_database.txt 
> ${working_directory_path}/project_database.txt" 
echo ${cmd_cat}
eval ${cmd_cat}

gunzip ${Name_of_project}_R1.fq.gz
gunzip ${Name_of_project}_R2.fq.gz

#run swarm with the different D values
#values passed by config file
#values="1 2 3 4"
for v in ${values}
do
	#mkdir cluster_d${v}
	#mkdir novel_d${v}
	# run the clustering on the REAL dataset
	echo "9) running swarm at ${v} thresholds"
	swarm_command="swarm -t ${num_threads} --append-abundance 1 -d ${v} 
	-o ${Name_of_project}_reads_cluseterd_with_phy_ITS_Swarmd${v}_coded.out 
	${working_directory_path}/temp_reads_plus_database.fasta" 
	echo ${swarm_command}
	eval ${swarm_command}
	wait
	
	#python rewrite the temp cluster names with the original old names
	# python rename the files so the names work in Swarm
	echo "10) python rename decode the cluster outfile"
	cmd_recode="python ${repository_path}/Python_ITS_scripts/parse_clusters.py 
	-i ${Name_of_project}_reads_cluseterd_with_phy_ITS_Swarmd${v}_coded.out 
	-d ${working_directory_path}/project_database.txt 
	-o ${working_directory_path}/${Name_of_project}_results/${Name_of_project}_reads_cluseterd_with_phy_ITS_Swarmd${v}.RESULTS" 
	echo ${cmd_recode}
	eval ${cmd_recode}
	wait
	

	# python graphs the  clussters
	echo "12) python to graphically represent clusters"
	cmd_summarise_graphs="python ${repository_path}/Python_ITS_scripts/draw_bar_chart_of_clusters.py 
	-i ${working_directory_path}/${Name_of_project}_results/${Name_of_project}_reads_cluseterd_with_phy_ITS_Swarmd${v}.RESULTS 
	-o temp" 
	echo ${cmd_summarise_graphs}
	eval ${cmd_summarise_graphs}
	wait
	
	# python find Phy species identified
	#echo "13) python to find Phy species identified "
	#cmd_what_phy_species="python ${repository_path}/Python_ITS_scripts/get_results_from_cluster_and_novel_clusterings.py --Name_of_project ${Name_of_project} --all_fasta ${working_directory_path}/temp_reads_plus_database_full_names.fasta -f ./${Name_of_project}_outfiles/fastqjoin.join.fasta --read_prefix ${read_prefix} -i ${working_directory_path}/${Name_of_project}_results/${Name_of_project}_reads_cluseterd_with_phy_ITS_Swarmd${v}.RESULTS 	-o ${working_directory_path}/${Name_of_project}_results/${Name_of_project}_Phytophthora_species_identified_swarm_${v}.txt --left ${Name_of_project}_R1.fq.gz --right ${Name_of_project}_R2.fq.gz" 
	#echo ${cmd_what_phy_species}
	#eval ${cmd_what_phy_species}
	#wait
	# python find Phy species identified
	echo "13) python to find Phy species identified "
	cmd_what_phy_species="python ${repository_path}/Python_ITS_scripts/get_results_from_cluster_and_novel_clusterings.py 
	--Name_of_project ${Name_of_project} 
	--all_fasta ${working_directory_path}/temp_reads_plus_database_full_names.fasta -f ./${Name_of_project}_outfiles/fastqjoin.join.fasta 
	--read_prefix ${read_prefix} 
	-i ${working_directory_path}/${Name_of_project}_results/${Name_of_project}_reads_cluseterd_with_phy_ITS_Swarmd${v}.RESULTS 
	-o ${working_directory_path}/${Name_of_project}_results/${Name_of_project}_Phytophthora_species_identified_swarm_${v}.txt 
	-v ${v}
	--threads ${num_threads}
	--left ${working_directory_path}/${Name_of_project}_R1.fq --right ${working_directory_path}/${Name_of_project}_R2.fq" 
	echo ${cmd_what_phy_species}
	eval ${cmd_what_phy_species}
done



gzip ${Name_of_project}_R1.fq
gzip ${Name_of_project}_R2.fq


rm *coded.out
rm temp*
rm project_database.txt



