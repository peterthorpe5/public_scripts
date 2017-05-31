#!/bin/bash
#$ -cwd
#Abort on any error,
#set -e

#echo Running on $HOSTNAME
#echo Current PATH is $PATH
#source $HOME/.bash_profile
temp="
Phytophthora_infestans_braker1.9_20170403
Phytophthora_infestans_RNAseq_short_train20170321
Phytophthora_infestans_RNAseq_short_train20170323
Phytophthora_infestans_RNAseq_train20170321"

################################################################
# Variables: FILLL IN DOWN TO THE END OF VARIABLES
Phy_dir=$HOME/scratch/tree_health/ITS_ratio/WORKED_Phytophthora_infestans.ASM14294v1.31
# known_fa="${Phy_dir}/P_infestans_genes_for_Pete_Thorpe.fasta"
known_fa="${Phy_dir}/Pi_T30_4.AA.fasta"
known_fa_nucl="${Phy_dir}/Pi_T30_4nt.fa"

Aug_species="Phytophthora_infestans_RNAseq_train20170321"

prefix="Pinf"
# default name
test_fa="aa.fa"
min_len_gene="20"
threads=8
python_directory=$HOME/public_scripts/gene_model_testing
Working_directory=${Phy_dir}/repeat_masking/${Aug_species}_repeat
test_gff="${Working_directory}/${Aug_species}.gff"
# for the repeat masking and GFF I used a altered gene name version
genome="${Phy_dir}/repeat_masking/Pi_alt.fasta"
#genome="${Phy_dir}/Phytophthora_infestans.ASM14294v1.31.fa"

# FOR HGT
# tax_filter_out is the phylum your beast lives in, or clade if you want to get a more refined HGT result
# tax_filter_up_to e.g. metazoan = tax_filter_up_to
# for aphid: #tax_filter_out=6656 #tax_filter_up_to=33208 
# for nematodes, #tax_filter_out=6231 ,#tax_filter_up_to=33208 

# for Phytophthora
T_30_4=403677
species_tx_id=4787
tax_filter_out=4783
#  Stramenopiles - heterokonts
tax_filter_up_to=33634 

# If you want to run transrate to get the RNAseq read mapping to gene 
# fill these out. Else, just let it fail, as it is the last step.

left_reads="${Phy_dir}/RNAseq_reads/R1.fq.gz"
right_reads="${Phy_dir}/RNAseq_reads/R2.fq.gz"

# END OF USER VARIABLES. NOTHING TO FILL IN FROM HERE.
#######################################################################
export Phy_dir
export known_fa
export known_fa_nucl
export prefix
export test_fa
export min_len_gene
export threads
export python_directory
export Working_directory
export test_gff
export genome
export T_30_4
export species_tx_id
export tax_filter_out
export tax_filter_up_to
export left_reads
export right_reads

#rm -rf ${Working_directory}
#mkdir ${Working_directory}
cd ${Working_directory}




# run Augustus with the desired models. 
echo "running Augustus"
aug="$HOME/scratch/Downloads/augustus-3.2.1/bin/augustus 
	--singlestrand=true 
	--genemodel=complete 
	--protein=on 
	--start=on 
	--alternatives-from-evidence=true 
	--stop=on 
	--cds=on 
	--introns=on 
	--gff3=off 
	--extrinsicCfgFile=/home/pt40963/Downloads/augustus-3.0.3/config/species/Myzus_cerasi/extrinsic.Myzus_cerasi.cfg 
	--hintsfile=${Phy_dir}/RNAseq_reads/hints_RNAseq.gff 
	--outfile=${test_gff} 
	--noInFrameStop=true 
	--species=${Aug_species}
	${genome}"
echo ${aug}
eval ${aug}
wait

echo "Augustus done"

/${python_directory}/Gene_model_testing_Master.sh


