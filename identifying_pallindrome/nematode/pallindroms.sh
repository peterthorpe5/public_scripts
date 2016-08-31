#!/bin/bash
set -e
# nemtaode mitochondria pallindrom finding
cd $HOME/nematode/Mitochondria_20160628

fasta_filenames="
Ge_mtDNAI.fa
Gp_mtDNAI.fa
Gp_mtDNAIII.fa
Gp_mtDNAV.fa
Gr_mtDNAI.fa
Gr_mtDNAIII.fa
Gr_mtDNAV.fa
Ge_mtDNAII.fa
Gp_mtDNAII.fa
Gp_mtDNAIV.fa
Gp_mtDNAVI.fa
Gr_mtDNAII.fa
Gr_mtDNAIV.fa"


for f in ${fasta_filenames}
do
	echo "Running  ${f}"

	#cd ${new_name}
	#step1
	cmd="$HOME/misc_python/identifying_pallindrome_final/identify_pallidromes.py --mask_NNN True --max_pallindrome_len 55 --fasta ${f} --base_name ${f}" 
	echo ${cmd}
	eval ${cmd}
	wait
	#step2
	cmd2="$HOME/misc_python/identifying_pallindrome_final/pallindrome_repeat_reducer.py --base_name ${f}" 
	echo ${cmd2}
	eval ${cmd2}
	wait
	#step3
	#cmd3="$HOME/misc_python/identifying_pallindrome_final/genic_or_non-genic.py --base_name ${f} --gff $HOME/nematode/mitochondria/genbank_files/${f}.gff" 
	#echo ${cmd3}
	#eval ${cmd3}
	wait
	#step4
	cmd4="$HOME/misc_python/identifying_pallindrome_final/pallindrom_clustering_assessment.py --base_name ${f}" 
	echo ${cmd4}
	eval ${cmd4}
	cd random
	new_name=${f%.fa}
	mkdir outfiles_${new_name}
	mv ${f}* ./outfiles_${new_name}
	mv random_${f}* ./outfiles_${new_name}
	cd ..
	
done

