READme Fix_five_prime_CDS
=========================

basic usage:

``./Fix_five_prime_CDS.py`` -h 

Usage: Use as follows:

``./python Fix_five_prime_CDS.py`` -t trnascriptome --cds nt_coding_seq --bam index_sorted_bam_file.bam
    --gff if_you_have_one --exp outfile_name_for_expression_values (default: overall_reads_mapped_per_sequences.txt)
    --prot (protein_cds_not_currently_used) -o outfile

Fix five prime cds based on read depth coverage v0.1.0:
why? Sometime the 5 prime end is not correctly predicted. This is really important to our
research (signal peptides etc ... Can we imporve this?

Requires:
samtools
Biopython
numpy


steps:

1) Index your transcriptome and map your reads back to your transcriptome:
    bowtie-build -f test_sequences.fasta test_sequences
    bowtie -S -p 8 test_sequences -1 R1.fq -2 R2.fq test_mapped.sam
2) Samtool index your transcriptome
    samtools faidx transcriptome.fasta
3) convert sam to bam
    samtools view -S -b -o unsorted.bam test_mapped.sam
4) sort your bam out!
    samtools sort unsorted.bam sorted.bam
5) Index your sorted.bam
    samtools index sorted.bam

How?

This script looks at the number of reads that map per base for your transcript of interest.
If the number of mapped reads is greater than "threshold (default=10)" then it performs statistical
analysis on this to determine if the starting methionoine (as found in the predicted cds)
has significantly lower expression that the rest of the transcript.

If it has, then this "ATG" may not be the coreect starting postion. It then looks for the next "ATG"
and performs statistical analysis to determine if this is a sensible starting postiton.

TO DO:

look upstream in transcript if CDS does not start with ATG, based on expression, could there be a good
candidate?

FIX five prime start in genomic regions based on this methodology



Options:
  -h, --help            show this help message and exit
  -t FILE, --transcriptome=FILE
                        the transcriptome assembly .fasta
  --cds=FILE            the predicted cds from the transcriptome assembly
                        .fasta
  --gff=FILE            the predicted coordinate for the cds predictions .gff
  --prot=FILE           the predicted amino acid cds  .pep
  -g FILE, --genome=FILE
                        the genome sequence. Not currently used. TO DO
  --bam=FILE            the sorted, indexed bam file as a result of the reads
                        beingmapped back to the transcriptome  .bam
  --min_read_count=MIN_READ_COUNT
                        the min_read_count that a transcript must have before
                        it isconsidered for statistical analysis. Default =
                        100
  --min_max_cov_per_base=MIN_MAX_COV_PER_BASE
                        the min_max_cov_per_base that a transcript must have
                        before it isconsidered for statistical analysis.
                        Default = 30
  --std=STAND_DEV_THRESHOLD
                        If the expression of the start of the cds is
                        stand_dev_threshold +/- the mean the flag as to be
                        looked at. Default = 3
  --help_full=HELP_FULL
                        prints out a full description of this program
  --exp=FILE            Output filename for the overall number of readsthat
                        map to the sequence of interest. Default:
                        overall_reads_mapped_per_sequences.txt
  -o FILE, --out=FILE   Output filename (default: results.out)


