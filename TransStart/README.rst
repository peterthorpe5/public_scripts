TranStart
=========

Transcription start finder based on read depth coverage v0.1.0
THIS IS NOT finding the start of coding squence. but where transcription binds and starts

``python TransStart.py -g genome.fasta --bam index_sorted_bam_file.bam --walk 3 --interation_value 1 --gff genes.gff -o outfile.txt``


Requirements:
    you must have samtools in your PATH,	python 3, Biopython, numpy

steps:
======

1) Index your transcriptome and map your reads back to your genome:
    Use STAR: a spice aware aligner
2) STAR --runMode genomeGenerate --runThreadN 10 --limitGenomeGenerateRAM
    74554136874 --genomeDir /PATH_TO/M.cerasi/star_indicies
    --genomeFastaFiles /PATH_TO/M.cerasi/Mc_v1.fasta
3) STAR --genomeDir star_indicies/ --runThreadN 12
    --outFilterMultimapNmax 5 --outSAMtype BAM SortedByCoordinate
    --outFilterMismatchNmax 7 --readFilesCommand zcat
    --outFileNamePrefix Mc --readFilesIn /PATH_TO/M.cerasi/RNAseq/Mc_R1.fq.gz
    /PATH_TO/M.cerasi/RNAseq/Mc_R2.fq.gz
2) Index your genome:
4) sort your bam out!
    samtools -@ 12 sort unsorted.bam sorted.bam
5) Index your sorted.bam
    samtools index sorted.bam

How?
This script looks at the number of reads that map per base for your gene
of interest.
If the number of mapped reads is greater than "threshold", reads
mapped per current base, then it performs statistical
analysis on this to determine if the start of transcription


Options:
  -h, --help            show this help message and exit
  --gff=FILE            the predicted coordinate for the cds predictions .gff
  -g FILE, --genome=FILE
                        the genome sequence. Not currently used. TO DO
  -b FILE, --bam=FILE   the sorted, indexed bam file as a result of the reads
                        being mapped back to the transcriptome  .bam
  -s STAND_DEV_THRESHOLD, --std=STAND_DEV_THRESHOLD
                        If the expression of the current region is less then
                        this number of std away from the first exon mean.
                        Default = 4
  -w WALK, --walk=WALK  the number of bases to walk away to find expression
                        Default = 3
  --min_value=MIN_VALUE
                        the min expression to return the coordinates for
                        Default = 2
  -i INTERATION_VALUE, --interation_value=INTERATION_VALUE
                        size of window to walk Default = 1
  --help_full=HELP_FULL
                        prints out a full description of this program
  --keep_gene_depth=KEEP_GENE_DEPTH
                        keep the output of the depth for the genes
  --logger=FILE         Output logger filename. Default: outfile_std.log
  -o FILE, --out=FILE   Output filename (default: results.out)

  
 convert the GFF output to fasta:
 ================================
``python GFF_to_fasta.py --gff TransStart_walk3based_on_min_value.gff -g nGr.v1.1.fa -o GROS_TSS_predition_maxLen1000.fasta -x 1000 -i 20``

Warning: This is made for a specific purpose. With specific data
      If you are using it, check it works for you!!
Usage: Use as follows:

python GFF_to_fasts.py --gff transtart.gff -m min len -x max len
        -g genome.fasta -o UTR.fasta

Options:
  -h, --help            show this help message and exit
  --gff=FILE            the predicted coordinate for the cds predictions .gff
  -g FILE, --genome=FILE
                        the genome sequence.
  -m FILE, --min_length=FILE
                        min_length of seq to return
  -x FILE, --max_length=FILE
                        max_length of seq to return
  -i FILE, --into_TSS=FILE
                        into_TSS of seq to return to the binding theory
  -u FILE, --upstream=FILE
                        upstream of TSS to return
  -o FILE, --out=FILE   Output filename (transtart.UTR.fasta)
