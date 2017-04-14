#!/bin/bash
#Abort on any error,
set -e

#echo Running on $HOSTNAME
#echo Current PATH is $PATH
#source ~/.bash_profile

################################################################
# Variables: FILLL IN DOWN TO THE END OF VARIABLES

Phy_dir=$HOME/scratch/tree_health/ITS_ratio/WORKED_Phytophthora_infestans.ASM14294v1.31

known_fa="${Phy_dir}/tests.AA.fasta"
known_fa_nucl="${Phy_dir}/Pi_T30_4nt.fa"
# default name
test_fa="aa.fa"
threads=8
python_directory=$HOME/public_scripts/gene_model_testing
Working_directory=$HOME/scratch/Pi/
test_gff="${Phy_dir}/repeat_masking/Pi.models_RNAseq.v1.VERY.short_NOT_BOTH_STRANDS.gff"
# for the repeat masking and GFF I used a altered gene name version
genome="${Phy_dir}/repeat_masking/Pi_alt.fasta"
#genome="${Phy_dir}/Phytophthora_infestans.ASM14294v1.31.fa"

# FOR HGT
# tax_filter_out is the phylum your beast lives in, or clade if you want to get a more refined HGT result
# tax_filter_up_to e.g. metazoan = tax_filter_up_to
# for aphid: #tax_filter_out=6656 #tax_filter_up_to=33208 
# for nematodes, #tax_filter_out=6231 ,#tax_filter_up_to=33208 

# for Phytophthora
species_tx_id=4787
tax_filter_out=4783 
tax_filter_up_to=33208 

# If you want to run transrate to get the RNAseq read mapping to gene 
# fill these out. Else, just let it fail, as it is the last step.

left_reads="${Phy_dir}//RNAseq_reads/R1.fq.gz"
right_reads="${Phy_dir}//RNAseq_reads/R2.fq.gz"

# END OF VARIABLES !!!!!!!!!!!!!!!!!!!!!!
########################################################

cd ${Working_directory}

echo "For this program you need:
Genome tools. Blast. gffread (cufflinks). Python. Biopython. diamond"

# genome tools to check and reformat the gff - essential step!
echo "prepare the amino acid seq from the GFF. First check gff"
gt_cmd="gt gff3 -sort -tidy -addids -sortnum -fixregionboundaries 
	   -addintrons -force -o ${test_gff}_reformatted.gff3 ${test_gff}"
echo ${gt_cmd}
eval ${gt_cmd}
wait

# GT to get stats on the predicted models. 

gt="gt gff3 -sort -tidy -addintrons ${test_gff}_reformatted.gff3 | 
   gt stat -force -genelengthdistri -o augustus_genelengthdistri.STAT 
   > temp"
echo ${gt}
eval ${gt}
wait
gt="gt gff3 -sort -tidy -addintrons ${test_gff}_reformatted.gff3 | 
gt stat -force -genescoredistri -o augustus_genescoredistri.STAT 
> temp"
echo ${gt}
eval ${gt}
wait
gt="gt gff3 -sort -tidy -addintrons ${test_gff}_reformatted.gff3 
| gt stat -force -exonlengthdistri -o augustus_exonlengthdistri.STAT 
> temp"
echo ${gt}
eval ${gt}
wait
gt="gt gff3 -sort -tidy -addintrons ${test_gff}_reformatted.gff3 
| gt stat -force -exonnumberdistri -o augustus_exonnumberdistri.STAT 
> temp"
echo ${gt}
eval ${gt}
wait
gt="gt gff3 -sort -tidy -addintrons ${test_gff}_reformatted.gff3 
| gt stat -force -intronlengthdistri -o augustus_intronlengthdistri.STAT 
> temp"
echo ${gt}
eval ${gt}
wait
gt="gt gff3 -sort -tidy -addintrons ${test_gff}_reformatted.gff3 
| gt stat -force -cdslengthdistri -o augustus_cdslengthdistri.STAT 
> temp"
echo ${gt}
eval ${gt}
wait

mkdir gff_stats
mv *.STAT ./gff_stats

rm temp

# gffread to the the cds from the reformat the gff - essential step!
echo "getting the AA and nt cds from the genome"
gff_to_cds="gffread ${test_gff}_reformatted.gff3 
		   -g ${genome} -x nt.fa 
		   -y aa.fa"
echo ${gff_to_cds}
eval ${gff_to_cds}
wait

#make blastdb
echo "step1: make blastdb"
mkdb="makeblastdb -in ${test_fa} -dbtype prot"
echo ${mkdb}
eval ${mkdb}
wait

#blast to xml
echo "step2: blast to xml"
bl_p="blastp -db ${test_fa} -query ${known_fa} -evalue 0.00001 
	  -seg no -num_threads ${threads} 
      -outfmt 5 -out test_fa_vs_known_fa.xml"
echo ${bl_p}
eval ${bl_p}
wait

mkdir known_fa_all_hits

#blast 2 tab
echo "step2b: blast to tab - top hit only"
bl_p2="blastp -db ${test_fa} -query ${known_fa} -evalue 1e-10 
	  -seg no -max_target_seqs 1 -num_threads ${threads} 
	  -outfmt 7 -out test_fa_vs_known_fa.tab"
echo ${bl_p2}
eval ${bl_p2}
wait

# convert the xml
echo "step3: convert the xml file"
parse_xml="python ${python_directory}/BLAST_parser_return_hits_NAME_only.py 
		  -i test_fa_vs_known_fa.xml 
		  -o test_fa_vs_known_fa.xml.condensed.out"
echo ${parse_xml}
eval ${parse_xml}
wait

# get the matches seq for the tab blast to the effectors
echo "step4: get the seqs of the top hit."
get_seq="python ${python_directory}/Get_sequence_from_tab_blast.py 
	    -b test_fa_vs_known_fa.tab 
		--known_fa ${known_fa}
		-p ${test_fa} 
		-n ${test_fa}"
echo ${get_seq}
eval ${get_seq}
wait

# change to where the file have been put
cd $HOME/${Working_directory}/known_fa

filenames=*.fasta

for f in ${filenames}
do
	echo "Running muscle ${f}"
	cmd="$HOME/Downloads/muscle3.8.31_i86linux64 -in ${f} 
		-out ${f}_aligned.fasta -maxiters 5000 -maxtrees 15" 
	echo ${cmd}
	eval ${cmd}
	wait
done

filenames2=*_aligned.fasta
for file in ${filenames2}
do
	echo "Running muscle ${f}"
	cmd="$HOME/Downloads/muscle3.8.31_i86linux64 
		-in ${file} -out ${file}_refine.fasta 
		-refine" 
	echo ${cmd}
	eval ${cmd}
	wait
done
	
rm *_aligned.fasta

mkdir alignments
mv *_refine.fasta ./alignments

###########################################################################
# diamond blast aginst NR?
echo "running diamond-BLAST against NR"
diam_p="diamond blastp -p 16 --sensitive -e 0.00001 
	   -v -q aa.fa 
	   -d /mnt/shared/scratch/pt40963/blast_databases/nr.dmnd 
	   -a aa.fasta_vs_nr.da"
echo ${diam_p}
eval ${diam_p}
wait

echo "converting diamond-BLAST output"
diam_v="diamond view -a aa.fasta*.daa -f tab -o aa.fasta_vs_nr.tab"
echo ${diam_v}
eval ${diam_v}
wait

echo "adding tx_id and descriptions to diamond-BLAST output"
tax="python ~/misc_python/diamond_blast_to_kingdom/Diamond_blast_to_taxid_add_kingdom_add_species_description.py 
	-i aa.fasta_vs_nr.tab 
	-p /home/pt40963/Downloads/blast_databases
	-o aa.fasta_vs_nr_tax.tab"
echo ${tax}
eval ${tax}
wait

echo "predicting HGT"
HGT="python /home/pt40963/misc_python/Lateral_gene_transfer_prediction_tool/Lateral_gene_transfer_predictor.py 
		-i *_vs_nr_tax.tab 
		--tax_filter_out ${tax_filter_out} 
		--tax_filter_up_to ${tax_filter_up_to}
		-p /home/pt40963/Downloads/blast_databases -o LTG_results.out"

echo ${HGT}
eval ${HGT}
wait
#

#Filter taxomony commands:
echo "filtering blast results"
filter_top_blasts="python ~/misc_python/BLAST_output_parsing/top_BLAST_hit_filter_out_tax_id.py 
				  -i *_vs_nr_tax.tab 
				  -t ${tax_filter_out} 
				  -p /home/pt40963/Downloads/blast_databases 
				  -o top_not_phylum_${tax_filter_out}.hits"
echo ${filter_top_blasts}
eval ${filter_top_blasts}
wait

filter_species="python ~/misc_python/BLAST_output_parsing/top_BLAST_hit_filter_out_tax_id.py 
			   -i *_vs_nr_tax.tab 
			   -t ${species_tx_id}
			   -p /home/pt40963/Downloads/blast_databases 
			   -o top_not_species_tx_id_${species_tx_id}.hits"
echo ${filter_species}
eval ${filter_species}
wait


#################################################################################
# change back to the working directory
# Run transrate?
cd ${Working_directory}
echo "Running transrate to get RNAseq mapping and it will tell you 
statistically what genes may be fusions. "
tran="transrate --assembly nt.fa 
	 --left ${left_reads} 
	 --right ${right_reads}"
echo ${tran}
eval ${tran}

echo ................................................................ im done
