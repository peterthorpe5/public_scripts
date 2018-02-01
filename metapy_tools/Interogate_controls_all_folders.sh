
cd /home/pt40963/scratch/tree_health/Illumina_results_2017Nov_first_run/controls

# this is a list of folder with control results in
folders="G*1_RESULTS"

# iterate through all the control folder and run the script to get the number
# of mitchase against the control sequences.
for folder in ${folders}
do
	cd ./${folder}
	echo "NOW RUNNING:  ${folder}"
	cd ./*_PEAR
	pwd
	fasta=*.assembled.fastqdrep.vsearch.fasta
	makeblastdb -in ${fasta} -dbtype nucl
	wait
	blastn -db ${fasta} -max_target_seqs 100000 -num_threads 8 -query /home/pt40963/scratch/tree_health/Illumina_results_2017Nov_first_run/controls.fasta -evalue 10 -outfmt 5 -out assembled_vs_control.xml 
	wait
	# creat a mismatch alignment and return coordinates of the mismatches
	python ~/public_scripts/Sanger_read_metagenetics/bin/BLAST_xml_parser_mismatch_coordinate.py -i assembled_vs_control.xml  -e 1e-5 -m 12 -o ${folder}_MISMATCHES.fastqdrep.vsearch.txt
	#value="0 1 2 3"
	# 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33"
	value="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
	for v in ${value}
	do
		echo "filtering:  ${v}"
		cmd="python ~/public_scripts/Sanger_read_metagenetics/bin/BLAST_xml_parser_for_MCL_exact.py 
		-i assembled_vs_control.xml
		-m ${v} 
		-o assembled_vs_control.exact.${v}_mismatch.txt"
		echo ${cmd}
		eval ${cmd}
		wait
		echo "c1  ${v}  mismatches:">> Number_found_using_dereplicated_reads_c1.txt
		cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "Control_05" >> Number_found_using_dereplicated_reads_c1.txt
		echo "c2  ${v}  mismatches::">> Number_found_using_dereplicated_reads_c2.txt
		cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "Control_14" >> Number_found_using_dereplicated_reads_c2.txt
		echo "c3  ${v}  mismatches::">> Number_found_using_dereplicated_reads_c3.txt
		cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "Control_22" >> Number_found_using_dereplicated_reads_c3.txt
		echo "c4  ${v}  mismatches::">> Number_found_using_dereplicated_reads_c4.txt
		cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "Control_42" >> Number_found_using_dereplicated_reads_c4.txt
		wait
		cmd2="python ~/public_scripts/Sanger_read_metagenetics/bin/BLAST_xml_parser_for_MCL.py 
		-i assembled_vs_control.xml 
		-m ${v} 
		-o assembled_vs_control.${v}_mismatch_or_less.txt"
		echo ${cmd2}
		eval ${cmd2}
	done
	rm -rf assembled_vs_control.xml
	find . -size 0 -delete
	rm -rf results_with_dereplicated_reads_vsearch
	mkdir results_with_dereplicated_reads_vsearch
	
	mv *.txt ./results_with_dereplicated_reads_vsearch
	#####################################################################################################
	
	fasta=*.assembled.fastqfor_swarm.fasta
	makeblastdb -in ${fasta} -dbtype nucl
	wait
	blastn -db ${fasta} -max_target_seqs 100000 -num_threads 8 -query /home/pt40963/scratch/tree_health/Illumina_results_2017Nov_first_run/controls.fasta -evalue 10 -outfmt 5 -out assembled_vs_control.xml 
	wait
	python ~/public_scripts/Sanger_read_metagenetics/bin/BLAST_xml_parser_mismatch_coordinate.py -i assembled_vs_control.xml  -e 1e-5 -m 12 -o ${folder}_MISMATCHES.fastqfor_swarm.txt
	for v in ${value}
	do
		echo "filtering:  ${v}"
		cmd="python ~/public_scripts/Sanger_read_metagenetics/bin/BLAST_xml_parser_for_MCL_exact.py 
		-i assembled_vs_control.xml
		-m ${v} 
		-o assembled_vs_control.exact.${v}_mismatch.txt"
		echo ${cmd}
		eval ${cmd}
		wait
		echo "c1  ${v}  mismatches:">> Number_found_using_dereplicated_reads_c1.txt
		cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "Control_05" >> Number_found_using_dereplicated_reads_c1.txt
		echo "c2  ${v}  mismatches::">> Number_found_using_dereplicated_reads_c2.txt
		cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "Control_14" >> Number_found_using_dereplicated_reads_c2.txt
		echo "c3  ${v}  mismatches::">> Number_found_using_dereplicated_reads_c3.txt
		cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "Control_22" >> Number_found_using_dereplicated_reads_c3.txt
		echo "c4  ${v}  mismatches::">> Number_found_using_dereplicated_reads_c4.txt
		cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "Control_42" >> Number_found_using_dereplicated_reads_c4.txt
		wait
		cmd2="python ~/public_scripts/Sanger_read_metagenetics/bin/BLAST_xml_parser_for_MCL.py 
		-i assembled_vs_control.xml 
		-m ${v} 
		-o assembled_vs_control.${v}_mismatch_or_less.txt"
		echo ${cmd2}
		eval ${cmd2}
	done
	rm -rf assembled_vs_control.xml
	find . -size 0 -delete
	rm -rf results_with_dereplicated_reads_swarm
	mkdir results_with_dereplicated_reads_swarm
	
	mv *.txt ./results_with_dereplicated_reads_swarm	
	#####################################################################################################
	echo "now running with all reads!"
	
	fasta_full=*.assembled.fastq.bio.chopped.fasta
	
	makeblastdb -in ${fasta_full} -dbtype nucl
	
	blastn -db ${fasta_full} -max_target_seqs 100000 -num_threads 8 -query /home/pt40963/scratch/tree_health/Illumina_results_2017Nov_first_run/controls.fasta -evalue 10 -outfmt 5 -out assembled_full_vs_control.xml
	wait
	python ~/public_scripts/Sanger_read_metagenetics/bin/BLAST_xml_parser_mismatch_coordinate.py -i assembled_full_vs_control.xml -e 1e-5 -m 12 -o ${folder}_MISMATCHES.full.txt
	
	for v in ${value}
	do
		echo "filtering:  ${v}"
		echo "${fasta}"
		cmd="python ~/public_scripts/Sanger_read_metagenetics/bin/BLAST_xml_parser_for_MCL_exact.py 
		-i assembled_full_vs_control.xml
		-m ${v} 
		-o assembled_full_vs_control.exact.${v}_mismatch.txt"
		echo ${cmd}
		eval ${cmd}
		wait
		echo "c1  ${v}  mismatches:">> Number_found_using_all_reads_c1.txt
		cat assembled_full_vs_control.exact.${v}_mismatch.txt | grep -c "Control_05" >> Number_found_using_all_reads_c1.txt
		echo "c2  ${v}  mismatches::">> Number_found_using_all_reads_c2.txt
		cat assembled_full_vs_control.exact.${v}_mismatch.txt | grep -c "Control_14" >> Number_found_using_all_reads_c2.txt
		echo "c3  ${v}  mismatches::">> Number_found_using_all_reads_c3.txt
		cat assembled_full_vs_control.exact.${v}_mismatch.txt | grep -c "Control_22" >> Number_found_using_all_reads_c3.txt
		echo "c4  ${v}  mismatches::">> Number_found_using_all_reads_c4.txt
		cat assembled_full_vs_control.exact.${v}_mismatch.txt | grep -c "Control_42" >> Number_found_using_all_reads_c4.txt
		wait
		cmd2="python ~/public_scripts/Sanger_read_metagenetics/bin/BLAST_xml_parser_for_MCL.py 
		-i assembled_full_vs_control.xml 
		-m ${v} 
		-o assembled_full_vs_control.${v}_mismatch_or_less.txt"
		echo ${cmd2}
		eval ${cmd2}
		wait
	done
	
	find . -size 0 -delete
	rm -rf assembled_full_vs_control.xml
	rm -rf results_with_all_reads
	mkdir results_with_all_reads
	wait
	mv *.txt ./results_with_all_reads
	cd ../../
done
