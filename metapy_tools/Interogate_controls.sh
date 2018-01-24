
cd $HOME/scratch/tree_health/Illumina_results_2017Nov_first_run/version_0.001_database

folders="GC1-0x_S82_L001_RESULTS
GC1-100x_S11_L001_RESULTS
GC1-10x_S94_L001_RESULTS
GC2-0x_S23_L001_RESULTS
GC2-100x_S47_L001_RESULTS
GC2-10x_S35_L001_RESULTS
GC3-0x_S59_L001_RESULTS
GC3-100x_S83_L001_RESULTS
GC3-10x_S71_L001_RESULTS
GL4-0x_S95_L001_RESULTS
GL4-100x_S24_L001_RESULTS
GL4-10x_S12_L001_RESULTS
GL5-0x_S36_L001_RESULTS
GL5-100x_S60_L001_RESULTS
GL5-10x_S48_L001_RESULTS
GL6-0x_S72_L001_RESULTS
GL6-100x_S96_L001_RESULTS
GL6-10x_S84_L001_RESULTS"

cd ./GC1-0x_S82_L001_RESULTS

cd ./GC1-0x_S82_L001_PEAR

fasta=*.assembled.fastqdrep.vsearch.fasta

makeblastdb -in ${fasta} -dbtype nucl
 
blastn -db ${fasta} -num_threads 8 -query $HOME/scratch/tree_health/Illumina_results_2017Nov_first_run/controls.fasta -evalue 1e-5 -outfmt 5 -out assembled_vs_control.xml 
 
value="0 1 2"
# 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33"
#value="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33"


for v in ${value}
do
	echo "filtering:  ${v}"
	echo "${fasta}"
	cmd="python ~/public_scripts/Sanger_read_metagenetics/bin/BLAST_xml_parser_for_MCL_exact.py 
	-i assembled_vs_control.xml
	-m ${v} 
	-o assembled_vs_control.exact.${v}_mismatch.txt"
	echo ${cmd}
	eval ${cmd}
	wait
	echo "c1  ${v}  mismatches:">> Number_found_using_dereplicated_reads_c1.txt
	cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "C1" >> Number_found_using_dereplicated_reads_c1.txt
	echo "c2  ${v}  mismatches::">> Number_found_using_dereplicated_reads_c2.txt
	cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "C2" >> Number_found_using_dereplicated_reads_c2.txt
	echo "c3  ${v}  mismatches::">> Number_found_using_dereplicated_reads_c3.txt
	cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "C3" >> Number_found_using_dereplicated_reads_c3.txt
	echo "c4  ${v}  mismatches::">> Number_found_using_dereplicated_reads_c4.txt
	cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "C4" >> Number_found_using_dereplicated_reads_c4.txt
	echo "c5  ${v}  mismatches::">> Number_found_using_dereplicated_reads_c5.txt
	cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "C5" >> Number_found_using_dereplicated_reads_c5.txt
	wait
	cmd2="python ~/public_scripts/Sanger_read_metagenetics/bin/BLAST_xml_parser_for_MCL.py 
	-i assembled_vs_control.xml 
	-m ${v} 
	-o assembled_vs_control.${v}_mismatch_or_less.txt"
	echo ${cmd2}
	eval ${cmd2}
done

rm assembled_vs_control.xml
find . -size 0 -delete
rm -rf results_with_dereplicated_reads
mkdir results_with_dereplicated_reads

mv *.txt ./results_with_dereplicated_reads

#####################################################################################################
echo "now running with all reads!"

fasta=*.assembled.fastq.bio.chopped.fasta

makeblastdb -in ${fasta} -dbtype nucl
 
blastn -db ${fasta} -num_threads 8 -query $HOME/scratch/tree_health/Illumina_results_2017Nov_first_run/controls.fasta -evalue 1e-5 -outfmt 5 -out assembled_vs_control.xml 
 
for v in ${value}
do
	echo "filtering:  ${v}"
	echo "${fasta}"
	cmd="python ~/public_scripts/Sanger_read_metagenetics/bin/BLAST_xml_parser_for_MCL_exact.py 
	-i assembled_vs_control.xml
	-m ${v} 
	-o assembled_vs_control.exact.${v}_mismatch.txt"
	echo ${cmd}
	eval ${cmd}
	wait
	echo "c1  ${v}  mismatches:">> Number_found_using_all_reads_c1.txt
	cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "C1" >> Number_found_using_all_reads_c1.txt
	echo "c2  ${v}  mismatches::">> Number_found_using_all_reads_c2.txt
	cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "C2" >> Number_found_using_all_reads_c2.txt
	echo "c3  ${v}  mismatches::">> Number_found_using_all_reads_c3.txt
	cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "C3" >> Number_found_using_all_reads_c3.txt
	echo "c4  ${v}  mismatches::">> Number_found_using_all_reads_c4.txt
	cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "C4" >> Number_found_using_all_reads_c4.txt
	echo "c5  ${v}  mismatches::">> Number_found_using_all_reads_c5.txt
	cat assembled_vs_control.exact.${v}_mismatch.txt | grep -c "C5" >> Number_found_using_all_reads_c5.txt
	wait
	cmd2="python ~/public_scripts/Sanger_read_metagenetics/bin/BLAST_xml_parser_for_MCL.py 
	-i assembled_vs_control.xml 
	-m ${v} 
	-o assembled_vs_control.${v}_mismatch_or_less.txt"
	echo ${cmd2}
	eval ${cmd2}	
done

wait
find . -size 0 -delete
rm assembled_vs_control.xml
rm -rf results_with_all_reads
mkdir results_with_all_reads
wait
mv *.txt ./results_with_all_reads
