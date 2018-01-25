#!/bin/bash
#$ -cwd
#Abort on any error,
set -e

cd $HOME/public_scripts/generate_ITS1_database

# first save the Excel sheet as a text file. We can update this with PANDAS later
# Pandas, currently doesnt work for me. 

###############################################################################
# PUT the db file name here
DATA_TAB_FILE="Phytophthora_ITS1_DB_v0.003_20180111.txt"
version="0.003"
#
###############################################################################


# activate VM
echo "activating VM. You will have to change this path!"

source ~/misc_python/THAPBI/THAPBI-pycits/venv-THAPBI-pycits-py3/bin/activate


# make sure no dos line endings
dos2unix ${DATA_TAB_FILE}


DATA_TAB_FILE2="${DATA_TAB_FILE%.*}"
# folder for results to go into
rm -rf database_file
mkdir database_file

# STEP: 1
# download all the sequences. Be sure to check the warnings
echo "download all the sequences. Be sure to check the warnings"
download_cmd="python ./bin/generate_database_direct_download.py 
			-t ${DATA_TAB_FILE} 
			-o ./database_file/${DATA_TAB_FILE2}_full_len.fasta 
			> ./database_file/${DATA_TAB_FILE2}_full_len.WARNING"
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
		./database_file/${DATA_TAB_FILE2}_full_len.fasta 
		> ./database_file/check_me_for_completeITS_may_be_a_lie.txt"
echo ${HMM_cmd}
eval ${HMM_cmd}
wait

# STEP: 3a
# here hmmsearch looks for domian. If the whole seq matches the hmm models it
# wont put it into the tableout. So we fish for these first
# e.g. 
# >> 9_Phytophthora_irrigata_EU334634  
#   [No individual domains that satisfy reporting thresholds (although complete target did)]
echo "get perfect matches which dont get in the domain table output"
ITS_cmd="python ./bin/get_complete_matches_from_hmmsearch.py 
		-i ./database_file/${DATA_TAB_FILE2}_full_len.fasta	
		--HMM_std_file ./database_file/check_me_for_completeITS_may_be_a_lie.txt 
		-o ./database_file/check_me_for_completeITS_HITS_may_be_a_lie.fasta"
echo ${ITS_cmd}
eval ${ITS_cmd}
wait

# STEP: 3
echo "reqrite the database-fasta as ITS1 sequence only"
ITS_cmd="python ./bin/get_DOMAIN_region_I_want_from_fasta_Nucleotide.py 
		-i ./database_file/${DATA_TAB_FILE2}_full_len.fasta	
		--hmm ./database_file/ITS1_hmm_domain_table.out 
		-o ./database_file/${DATA_TAB_FILE2}_ITS_only.fasta"
echo ${ITS_cmd}
eval ${ITS_cmd}
wait

# STEP: 4
echo "muscle align the db .fasta to see any mistakes"
align_cmd="muscle 
			-in ./database_file/${DATA_TAB_FILE2}_ITS_only.fasta 
			-out ./database_file/${DATA_TAB_FILE2}_ITS_only.fasta.aln "
echo ${align_cmd}
eval ${align_cmd}
wait

# STEP: 4b
echo "muscle align refine the db .fasta "
align_cmdr="muscle 
			-in ./database_file/${DATA_TAB_FILE2}_ITS_only.fasta.aln 
			-out ./database_file/${DATA_TAB_FILE2}_ITS_only.fasta.aln.refine
			-refine"
echo ${align_cmdr}
eval ${align_cmdr}
wait

mv_files="mv  ./database_file/${DATA_TAB_FILE2}_ITS_only.fasta.aln.refine
	     ./database_file/${DATA_TAB_FILE2}_ITS_only.fasta.aln"
echo ${mv_files}
eval ${mv_files}
wait

# STEP: 5
echo "line wrap the fasta file"
ITS_cmd="python ./bin/rewrite_as_fasta.py 
		-i ./database_file/${DATA_TAB_FILE2}_ITS_only.fasta	
		-o ./database_file/Phytophora_ITS_database_v${version}.fasta"
echo ${ITS_cmd}
eval ${ITS_cmd}
wait

# STEP: 5b
echo "check for illegal characters (Only ATCG allowed"
ITS_cmd="python ./bin/Check_protein_seq_for_ILLEGAL_charcter.py 
		-i ./database_file/Phytophora_ITS_database_v${version}.fasta	
		-o temp.fasta >  ./database_file/Duplicate.WARNINGS"
echo ${ITS_cmd}
eval ${ITS_cmd}
wait

# STEP: 6
echo "Make sure no duplicates and illegal characters are still present"
ITS_cmd="python ./bin/Check_protein_seq_for_ILLEGAL_charcter.py 
		-i temp.fasta	
		-o temp_v${version}.fasta >  ./database_file/Duplicate2.WARNINGS"
echo ${ITS_cmd}
eval ${ITS_cmd}
wait

# STEP: 7
echo "line wrap the fasta file"
ITS_cmd="python ./bin/rewrite_as_fasta.py 
		-i temp_v${version}.fasta	
		-o Phytophora_ITS_database_v${version}.fasta"
echo ${ITS_cmd}
eval ${ITS_cmd}
wait

rm temp.fasta
rm temp_v${version}.fasta







