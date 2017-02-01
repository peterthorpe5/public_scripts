#!/bin/bash
# put out files in current directory
#$ -cwd
#$ -l hostname="n13*"
# Deliver notifications to the following address
# Send notifications when the job begins (b), ends (e), or is aborted (a)

#Abort on any error,
set -e
#(and thus quit the script right away)

echo Running on $HOSTNAME
echo Current PATH is $PATH
#source ~/.bash_profile
#export PATH=$HOME/bin:$PATH
#echo Revised PATH is $PATH

#which blastx

cd /home/pt40963/scratch/tree_health/nicotianae
export TMP=/home/pt40963/scratch/${USER}_${JOB_ID}
#prepare the data
# grep for gene field in the gff
#mv GCA_001482985.1_ASM148298v1_genomic.gff p.nicotiana_genome.gff
#cat p.nicotiana_genome.gff | grep "gene" | grep -v "mRNA" > p.nicotiana_genome.gene.gff
#mv GCA_001482985.1_ASM148298v1_genomic.fna p.nicotiana_genome.fasta

# blast to get representative ITS regions.
#I changed the names - badly!
#makeblastdb -in p.nicotiana_genome.fasta -dbtype nucl
#blastn -query P.infestnas_ITS.fasta -db p.nicotiana_genome.fasta -outfmt 6 -out n.Pi_ITS_vs_p.nictoiana.out

# prepare ITS gff. python:

python ~/misc_python/THAPBI/ITS_region_genomic_coverage/generate_ITS_GFF.py --blast n.Pi_ITS_vs_p.nictoiana.out --prefix P.nictoiana -o P.nictoiana.ITS.GFF

#reads 
 #wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR219/006/SRR2198696/SRR2198696_1.fastq.gz
 # wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR219/006/SRR2198696/SRR2198696_2.fastq.gz


 #quality trim the reads
#java -jar /home/pt40963/Downloads/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 12 -phred33 SRR2198696_1.fastq.gz SRR2198696_2.fastq.gz R1.fq.gz unpaired_R1.fq.gz R2.fq.gz R2_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 HEADCROP:9 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:51


mkdir $TMP

# index genome
bowtie2-build -f p.nicotiana_genome.fasta Pn

# randomly assign mutliple mapping reads ...  http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#bowtie2-options-p
#bowtie2 --very-sensitive --non-deterministic --seed 1 --no-mixed --no-unal -p 8 -x Pi -1 R1.fq.gz -2 R2.fq.gz -S P.nicotiana.sam

#pipe SAM output on stdout into samtools to make it into BAM
#TODO: Put the temp unsorted BAM file on the local hard drive scratch space under /mnt/scratch
bowtie2 --very-sensitive --non-deterministic --seed 1 --no-mixed --no-unal -p 16 -x Pn -1 R1.fq.gz -2 R2.fq.gz | samtools view -S -b -o $TMP/tmp_unsorted.bam -


# convert to sorted bam.
#samtools view -@ 8 -S -b -o tmp_unsorted.bam P.nicotiana.sam

samtools sort -@ 16 $TMP/tmp_unsorted.bam P.nicotiana
samtools index P.nicotiana.bam

# Clean up, doing this last in case fails
rm $TMP/tmp_unsorted.bam
# WARNING - Blindly deleting folders based on an encironemnt variable is DANGEROUS!
# There I am NOT doing: rm -rf $TMP
# rmdir will only work if the directory is emtpy (and it should be)
rmdir $TMP

# use bedtools to get the number of reads that map to specific regions

bedtools multicov -bams P.nicotiana.bam -bed p.nicotiana_genome.gene.gff > Phynico_t30_genomic.genes.cov

bedtools multicov -bams P.nicotiana.bam -bed p.nicotiana_genome.ITS_hits.gff > Phynico_t30_genomic.ITS.cov


# just get the values using cut

cat Phynico_t30_genomic.genes.cov | cut -f10 > Phynico_t30_genomic.genes.cov.values

cat Phynico_t30_genomic.ITS.cov | cut -f10 > Phynico_t30_genomic.ITS.cov.values


