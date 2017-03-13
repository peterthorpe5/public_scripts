#!/bin/bash
# put out files in current directory
#$ -cwd
#$ -l hostname="n13*"
# Deliver notifications to the following address
# Send notifications when the job begins (b), ends (e), or is aborted (a)
#dollar -M email.address@somewhere .... currently doesnt work
#dollar -m a 
#Abort on any error,
#set -euo pipefail
#(and thus quit the script right away)

echo Running on $HOSTNAME
echo Current PATH is $PATH
#source ~/.bash_profile

cd ~/scratch/tree_health/ITS_ratio
export TMP=~/scratch/${USER}_${JOB_ID}

##################################################################################################################################################################
# THESE VARIABLE NEED TO BE FILLED IN BY USER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

species=Phytophthora_kernoviae

genome_prefix=Phytophthora_kernoviae.GCA_000333075.1.31

genome_fasta=ftp://ftp.ensemblgenomes.org/pub/protists/release-31/fasta/phytophthora_kernoviae/dna/Phytophthora_kernoviae.GCA_000333075.1.31.dna.genome.fa.gz

genome_GFF=ftp://ftp.ensemblgenomes.org/pub/protists/release-31/gff3/phytophthora_kernoviae/Phytophthora_kernoviae.GCA_000333075.1.31.gff3.gz

read_1_link=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR639/SRR639379/SRR639379_1.fastq.gz

read_2_link=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR639/SRR639379/SRR639379_2.fastq.gz

read_1_b_link=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR278/008/SRR2785298/SRR2785298_1.fastq.gz
read_2_b_link=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR278/008/SRR2785298/SRR2785298_2.fastq.gz

trimmomatic_path=~/Downloads/Trimmomatic-0.32

SRA_prefix=SRR639379

repository_path=~/misc_python/THAPBI/ITS_region_genomic_coverage

num_threads=12


# NOTHING TO BE FILLED IN BY USER FROM HERE!!!!
##################################################################################################################################################################

export species
export genome_prefix
export genome_fasta
export genome_GFF
export read_1_link
export read_2_link
export trimmomatic_path
export SRA_prefix
export repository_path
export num_threads


#./ITS_genomic_coverage_ratio_finding.sh

cd Phytophthora_kernoviae.GCA_000333075.1.31

echo "MAPPING"
cmd_mapping="bowtie2 --very-sensitive --non-deterministic --seed 1 --no-mixed --no-unal -p ${num_threads} -x bowtie_index_files -1 R1.fq.gz -2 R2.fq.gz | samtools view -@ ${num_threads} -S -b -o $TMP/tmp_unsorted.bam -" 
echo ${cmd_mapping}
eval ${cmd_mapping}
echo "mapping fininshed"

# convert to sorted bam.
#samtools view -@ 8 -S -b -o tmp_unsorted.bam P.nicotiana.sam
echo "sort bam file"
cmd_sort="samtools sort -@ ${num_threads} $TMP/tmp_unsorted.bam ${genome_prefix}"
echo ${cmd_sort}
eval ${cmd_sort}

#index bam file
samtools index ${genome_prefix}.bam


# Clean up, doing this last in case fails
rm $TMP/tmp_unsorted.bam
# WARNING - Blindly deleting folders based on an encironemnt variable is DANGEROUS!
# There I am NOT doing: r m -r f $ T M P
# rmdir will only work if the directory is emtpy (and it should be)
rmdir $TMP

# get only the genes, not bothered about other stuff ...
echo "prepare GFF for genes only"
cmd_python_gene_to_gff="python ${repository_path}/get_genes_from_GFF.py --gff ${genome_prefix}.gff3 -o ${genome_prefix}.gene.gff"
echo ${cmd_python_gene_to_gff}
eval ${cmd_python_gene_to_gff}
 


# use bedtools to get the number of reads that map to specific regions

echo "bedtools count"
cmd_Busco_count="bedtools multicov -bams ${genome_prefix}.bam -bed ${genome_prefix}_BUSCO_GENES.gene.gff > ${genome_prefix}_BUSCO_GENES.gene.gff.genes.cov"
echo ${cmd_Busco_count}
eval ${cmd_Busco_count}

cmd_count="bedtools multicov -bams ${genome_prefix}.bam -bed ${genome_prefix}.gene.gff > ${genome_prefix}_genomic.genes.cov"
echo ${cmd_count}
eval ${cmd_count}

cmd_ITS_count="bedtools multicov -bams ${genome_prefix}.bam -bed ${genome_prefix}.ITS.consensus.GFF > ${genome_prefix}_genomic.ITS.cov"
echo ${cmd_ITS_count}
eval ${cmd_ITS_count}

echo counting done


# just get the values using cut
# remove RNA genes from the all gene values, as these are what we are measuring in the other
# dataset. Assuming they are annotated in the GFF. It not... results may be a bit skewed.


cut -f10 ${genome_prefix}_genomic.genes.cov > ${genome_prefix}_genomic.genes.cov.values
echo cut -f10 ${genome_prefix}_genomic.genes.cov to ${genome_prefix}_genomic.genes.cov.values

# for the ITS genes regions
cut -f10 ${genome_prefix}_genomic.ITS.cov >  ${genome_prefix}_genomic.ITS.cov.values
echo cut -f10 ${genome_prefix}_genomic.ITS.cov to  ${genome_prefix}_genomic.ITS.cov.values

#for the busco genes
cut -f10 ${genome_prefix}_BUSCO_GENES.gene.gff.genes.cov>  ${genome_prefix}_BUSCO_GENES.gene.gff.genes.cov.values
echo cut -f10 ${genome_prefix}_BUSCO_GENES.gene.gff.genes.cov to  ${genome_prefix}_BUSCO_GENES.gene.gff.genes.cov.values

# get stats summary of coverages
python ${repository_path}/summary_stats.py --ITS ${genome_prefix}_genomic.ITS.cov.values --GFF ${genome_prefix}.ITS.consensus.GFF --all_genes_cov ${genome_prefix}_genomic.genes.cov.values -o ${genome_prefix}_stats_all_genes_versus_ITS.out
echo python ${repository_path}/summary_stats.py --ITS ${genome_prefix}_genomic.ITS.cov.values --GFF ${genome_prefix}.ITS.consensus.GFF --all_genes_cov ${genome_prefix}_genomic.genes.cov.values -o ${genome_prefix}_stats_all_genes_versus_ITS.out

# get stats summary of coverages with BUSCO genes
python ${repository_path}/summary_stats.py --ITS ${genome_prefix}_genomic.ITS.cov.values --GFF ${genome_prefix}.ITS.consensus.GFF --all_genes_cov ${genome_prefix}_BUSCO_GENES.gene.gff.genes.cov.values -o ${genome_prefix}_stats_BUSCO_versus_ITS.out
echo python ${repository_path}/summary_stats.py --ITS ${genome_prefix}_genomic.ITS.cov.values --GFF ${genome_prefix}.ITS.consensus.GFF --all_genes_cov ${genome_prefix}_BUSCO_GENES.gene.gff.genes.cov.values -o ${genome_prefix}_stats_BUSCO_versus_ITS.out


