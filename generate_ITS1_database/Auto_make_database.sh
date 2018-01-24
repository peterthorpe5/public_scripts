#!/bin/bash
#$ -cwd
#Abort on any error,
set -e

cd $HOME/public_scripts/generate_ITS1_database


# first save the Excel sheet as a text file. We can update this with PANDAS later
# Pandas, currently doesnt work for me. 

# activate VM
echo "activating VM. You will have to change this path!"

source ~/misc_python/THAPBI/THAPBI-pycits/venv-THAPBI-pycits-py3/bin/activate

# make sure no dos line endings
dos2unix Phytophthora_ITS1_DB_v0.003_20180124.txt

# folder for results to go into
rm -rf database_file
mkdir database_file

# STEP: 1
# download all the sequences. Be sure to check the warnings
echo "download all the sequences. Be sure to check the warnings"
download_cmd="python ./bin/generate_database_direct_download.py 
			-t Phytophthora_ITS1_DB_v0.003_20180124.txt 
			-o ./database_file/Phytophthora_FULL_LEN_DB_v0.003_20180124.fasta 
			> ./database_file/Phytophthora_FULL_LEN_DB_v0.003_20180124.WARNING"
echo ${download_cmd}
eval ${download_cmd}
wait

# to make the hmm  profile. First align and refine the region of interest.
# Then make the hmm on the alingment:
# hmmbuild --dna Phytophora_ITS_database_v0.003.hmm Phytophora_ITS_database_v0.003.aln.refine

# STEP: 2
echo "HMM to identify ITS1 coordinates: REQUIRE hmmsearch"
HMM_cmd="hmmsearch 
		--domE 1e-22
		--domtblout ./database_file/ITS1_hmm_domain_table.out 
		./bin/Phytophora_ITS_20180104.hmm
		./database_file/Phytophthora_FULL_LEN_DB_v0.003_20180124.fasta"
echo ${HMM_cmd}
eval ${HMM_cmd}
wait

# STEP: 3
echo "reqrite the database-fasta as ITS1 sequence only"
ITS_cmd="python ./bin/get_DOMAIN_region_I_want_from_fasta_Nucleotide.py 
		-i ./database_file/Phytophthora_FULL_LEN_DB_v0.003_20180124.fasta	
		--hmm ./database_file/ITS1_hmm_domain_table.out 
		-o ./database_file/Phytophthora_ITS1_DB_v0.003_20180111_ITS_only.fasta"
echo ${ITS_cmd}
eval ${ITS_cmd}
wait

# STEP: 4
echo "muscle align the db .fasta to see any mistakes"
align_cmd="muscle 
			-in ./database_file/Phytophthora_ITS1_DB_v0.003_20180111_ITS_only.fasta 
			-out ./database_file/Phytophthora_ITS1_DB_v0.003_20180111_ITS_only.fasta.aln "
echo ${align_cmd}
eval ${align_cmd}
wait


# STEP: 5
echo "line wrap the fasta file"
ITS_cmd="python ./bin/rewrite_as_fasta.py 
		-i ./database_file/Phytophthora_ITS1_DB_v0.003_20180111_ITS_only.fasta	
		-o ./database_file/Phytophora_ITS_database_v0.003.fasta"
echo ${ITS_cmd}
eval ${ITS_cmd}
wait

# STEP: 5
echo "line wrap the fasta file"
ITS_cmd="python ./bin/Check_protein_seq_for_ILLEGAL_charcter.py 
		-i ./database_file/Phytophora_ITS_database_v0.003.fasta	
		-o temp.fasta >  ./database_file/Duplicate.WARNINGS"
echo ${ITS_cmd}
eval ${ITS_cmd}
wait

# STEP: 6
echo "line wrap the fasta file"
ITS_cmd="python ./bin/Check_protein_seq_for_ILLEGAL_charcter.py 
		-i temp.fasta	
		-o temp_v0.003.fasta >  ./database_file/Duplicate2.WARNINGS"
echo ${ITS_cmd}
eval ${ITS_cmd}
wait

# STEP: 7
echo "line wrap the fasta file"
ITS_cmd="python ./bin/rewrite_as_fasta.py 
		-i temp_v0.003.fasta	
		-o Phytophora_ITS_database_v0.003.fasta"
echo ${ITS_cmd}
eval ${ITS_cmd}
wait

rm temp.fasta
rm temp_v0.003.fasta







