#!/bin/bash
#$ -cwd

cd $HOME/scratch/Mpersicae/genotypeO/Fix_5_prime/test_data

# get the bam file
wget https://dl.dropboxusercontent.com/u/68811826/test_mapped100.sorted.bam

# transrate command - this wraps Snap which is really fast
#transrate --assembly test_sequences.fasta --left ~/scratch/Mp_O/Mp_O_RNAseq/Mp_O_R1.fq.gz --right ~/scratch/Mp_O/Mp_O_RNAseq/Mp_O_R2.fq.gz --threads 6
# index the fa file
samtools faidx test_sequences.fasta
# get the mapped reads only. 
#samtools view -@ 4 -S -b -F 4 unsorted.bam >  mapped.bam
samtools sort -@ 4 mapped.bam test_mapped100.sorted
samtools index test_mapped100.sorted.bam

# the bam file is too big for github, so subsample at 50%
samtools view -@ 4 -b -s 0.50 mapped.bam > mapped.0.5subsample.bam
samtools sort -@ 4 mapped.0.5subsample.bam mapped.0.5subsample.sorted
samtools index mapped.0.5subsample.sorted.bam

# the bam file is too big for github, so subsample at 75%
samtools view -@ 4 -b -s 0.25 mapped.bam > mapped.0.25subsample.bam
samtools sort -@ 4 mapped.0.25subsample.bam mapped.0.25subsample.sorted
samtools index mapped.0.25subsample.sorted.bam

# the bam file is too big for github, so subsample at 95%
samtools view -@ 4 -b -s 0.05 mapped.bam > mapped.0.05subsample.bam
samtools sort -@ 4 mapped.0.05subsample.bam mapped.0.05subsample.sorted
samtools index mapped.0.05subsample.sorted.bam


#bam read to identify chimeras. The probabilty that it is one gene - not more than one fused. 
#~/transrate-tools/src/bam-read name_offile.bam chimera_indetification_outfile.csv 0.95

# run a test:
#cd $HOME/Mpersicae/genotypeO/Fix_5_prime/test_data
# full data
python ~/misc_python/Fix_five_prime/Fix_five_prime_CDS.py -t test_sequences.fasta --bam test_mapped100.sorted.bam --cds cds.test_sequences.fasta --std 3 -o test.sd3.fasta > test.sd3.info
python ~/misc_python/Fix_five_prime/Fix_five_prime_CDS.py -t test_sequences.fasta --bam test_mapped100.sorted.bam --cds cds.test_sequences.fasta --std 5 -o test.sd5.fasta > test.sd5.info


python ~/misc_python/Fix_five_prime/Fix_five_prime_CDS.py -t test_sequences.fasta --bam mapped.0.5subsample.sorted.bam --cds cds.test_sequences.fasta --std 3 -o test.sub.sample0.5.sd3.fasta > test.sub.sample0.5.sd3.info
python ~/misc_python/Fix_five_prime/Fix_five_prime_CDS.py -t test_sequences.fasta --bam mapped.0.5subsample.sorted.bam --cds cds.test_sequences.fasta --std 5 -o test.sub.sample0.5.sd5.fasta > test.sub.sample0.5.sd5.info


python ~/misc_python/Fix_five_prime/Fix_five_prime_CDS.py -t test_sequences.fasta --bam mapped.0.25subsample.sorted.bam --cds cds.test_sequences.fasta --std 3 -o test.sub.sample0.25.sd3.fasta > test.sub.sample0.25.sd3.info
python ~/misc_python/Fix_five_prime/Fix_five_prime_CDS.py -t test_sequences.fasta --bam mapped.0.25subsample.sorted.bam --cds cds.test_sequences.fasta --std 5 -o test.sub.sample0.25.sd5.fasta > test.sub.sample0.25.sd5.info


python ~/misc_python/Fix_five_prime/Fix_five_prime_CDS.py -t test_sequences.fasta --bam mapped.0.05subsample.sorted.bam --cds cds.test_sequences.fasta --std 3 -o test.sub.sample0.05.sd3.fasta > test.sub.sample0.05.sd3.info
python ~/misc_python/Fix_five_prime/Fix_five_prime_CDS.py -t test_sequences.fasta --bam mapped.0.05subsample.sorted.bam --cds cds.test_sequences.fasta --std 5 -o test.sub.sample0.05.sd5.fasta > test.sub.sample0.05.sd5.info
