
please cite:
============

``Thorpe, P., Escudero-Martinez, C.M., Cock, P.J., Eves-van den Akker, S. and Bos, J.I., 2018. Shared transcriptional control and disparate gain and loss of aphid parasitism genes. Genome biology and evolution, 10(10), pp.2716-2733.``


READme info for script to Identify putative HGT events
======================================================

basic usage:

``./Lateral_gene_transfer_predictor.py`` -h 



script to open up a tab blast output and generate Alien index scores,
this finds the kingdom that blast mtach has been assigned.
script: Parse the NCBI taxonomic database node.dmp to get the
taxonomic relationship

author: Peter Thorpe September 2015. The James Hutton Insitute, Dundee, UK.

./Lateral_gene_transfer_predictor.py -h
Usage: Use as follows:

$ ``./Lateral_gene_transfer_predictor.py -i blast_w_tax_id.tab --tax_filter_out 6656 (e.g.arthropoda) --tax_filter_up_to 33208 (e.g. metazoan) -o LTG_results.out``


- for info: taxid - 6231 (nematoda)

Tax databse from NCBI is require. Download, unzip, and use -p /PATH/TO/   scripts will find them from here.

    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    
	tar -zxvf taxdump.tar.gz


e.g. 6656 = filter_out_tax_id --tax_filter_out
e.g. 33208 = Metazoa   -  for me this is the tax id I want to go up to --tax_filter_up_to

What:
To determine Lateral gene transfer event (LGT). An alien index score needs to be generated. Score > 45
is a candidate LGT - however, it could be contamination. The user will have to decide.


How:

Alien index: (plagerised from Seb Eves Van Den Akker et al, Naccubus paper....)
http://gbe.oxfordjournals.org/content/6/9/2181.short
taxonomic origin: Bacteria, Fungi, Plantae, Metazoa, and Other (including protists), and
the best expect value from each group was used to calculate an Alien Index (AI) as given
by the formula AI = log((Best E-value for Metazoa) + e-200) - log((Best E-value for Non-
Metazoa) + e-200). If no hits were found, the E-value was set to 1. Thus, AI is allowed to
vary in the interval between +460 and -460, being positive when top non-metazoan hits
yielded better E-values than the top metazoan ones. Entries with incomplete taxonomy,
such as X-ray structures, or hits belonging to the same phylum as the query sequence (
i.e Rotifera, filter_out_tax_id or Nematoda), were excluded from the analysis


    BLAST DATA should be formatted as:

1) qseqid = Query Seq-id (ID of your sequence)
2) sseqid = Subject Seq-id (ID of the database hit)
3) pident = Percentage of identical matches
4) length = Alignment length
5) mismatch = Number of mismatches
6) gapopen = Number of gap openings
7) qstart = Start of alignment in query
8) qend = End of alignment in query
9) sstart = Start of alignment in subject (database hit)
10) send = End of alignment in subject (database hit)
11) evalue = Expectation value (E-value)
12) bitscore = Bit score
13) salltitles = TOP description of the blast hit
14) staxids = tax_id
15) scientific_name
16) scomnames = common_name
17) sskingdoms = kingdom



MORE INFO:

Alien index:  (http://www.sciencemag.org/content/suppl/2008/05/29/320.5880.1210.DC1/Gladyshev.SOM.pdf)
    Massive Horizontal Gene Transfer in Bdelloid Rotifers :
    taxonomic origin: Bacteria, Fungi, Plantae, Metazoa, and Other (including protists), and
    the best expect value from each group was used to calculate an Alien Index (AI) as given
    by the formula AI = log((Best E-value for Metazoa) + e-200) - log((Best E-value for Non-
    Metazoa) + e-200). If no hits were found, the E-value was set to 1. Thus, AI is allowed to
    vary in the interval between +460 and -460,being positive when top non-metazoan hits
    yielded better E-values than the top metazoan ones. 

    Based on AI, genes were classified as foreign (AI>45), indeterminate (0<AI<45), or metazoan (AI<0)

There will be Four outfiles, using the -o prefix.

1) prefix name :
 This will contain the best metazoan and non metazoan hit per blast query. Depending on your filters specified.

2) prefix_precursor_value_temp :
 this file contain the precursor value for each of those if out_file(1). This is used for the next file

3) prefix_Alien_index.out :
 This file contains all the AI scores for each BLAST query. In one long line the best metazoanan
 and non-metazoan hits can be seen. A final comment is added if the programs believes there to be a HGT/LGT event,
 or if this is contamination.

 Contamination is based on a high AI score and greater than 70% identity to a non-metazoan. This can ofcourse be changed to suit the user.

4) Perfix_LGT_candifates.out :
 This file contains all scores greater than 0. The final comments box is a note to say if it think it is potential
 contamination or if it may be a HGT/LTG event.



Options:
  -h, --help            show this help message and exit
  -i FILE, --in=FILE    the tab output from blast/diamond. This must have
                        tax_id info in this!!If you use diamond, please get
                        this info using 'add_taxonomic_info_to_tab_output.py'
  -p PATH, --path=PATH  Directory containing relevant taxonomy/database files
                        Default is the current working directory. This is not
                        used with the main input and output filenames.
						
  --pi=PI               this is a threshold for determining likely
                        contanimants. e.g. if it is greater than pi percentage
                        identityt than it may be contanimantion.  or a very
                        recent HGT. Default = 70.
  --tax_filter_out=TAX_FILTER_OUT
                        The tax ID to filter out: for this analysis the Phylum
                        which your BEASTof interest if found. e.g. Aphids are
                        from Arthropoda, therefore this would be 6656, whihc
                        is the dwefault value. This will filter out all blast
                        hit which are from this phylum. It is possible to put
                        a species/kingdom tax_id in here ... whatever floats
                        your boat.
  --tax_filter_up_to=TAX_FILTER_UP_TO
                         The tax_id to 'walk up to', to determine assignment.
                        By default this is metazoa.The script work out the
                        best metazoan to non-metazoan hit. But this can be
                        altered if you wish to alter this
  --tax_coloumn=TAX_COLOUMN
                        the coloumn with the tax_id info. Defulat is 14(as
                        counted by a human/ not a computer
  -o FILE, --out=FILE   Output filename - default=
                        infile__tab_blast_LGT_result
					
					
Note: this script currently only ranges from -200 to +200. Not the range specified in their publication. 
Maybe an alterantive LOG is used.

TO DO:

This script does not yet bin the blast hits to kingdom. Im not entirely sure why this is done.


===================================================================================================

READme check contaminants_on_contigs
====================================

basic usage:

``./check_contaminants_on_contigs.py`` -h 


check_contaminants_on_contigs.py --gff ../augustus.gff3 -LTG LTG_LGT_candifates.out (default)

Title:
script to open gff and create a dictionary of {scaffold: set([gene1, gene2])
 this can then be used to see if all genes on a scaff are predicted to be HGT and therefore
 the whole scaffold is most likely contamination. 
 The script will output a file with contigs that only have contigs/scaffolds
 that are HGT/LTG genes

 
Tool to refine the HGT predicted gene based on RNAseq cov, genomic cov, exon number, percentage identity to best non-metazoan hit and AT content that differes from normal.

``python ~/misc_python/Lateral_gene_transfer_prediction_tool/check_contaminants_on_contigs.py`` --gff ../augustus.gff3 -LTG LTG_LGT_candifates.out (default)

``python ~/misc_python/Lateral_gene_transfer_prediction_tool/check_contaminants_on_contigs.py`` --bam sorted.bam --gff augustus.gff3 --LTG LTG_LGT_candifates_AI_30plus.out -s 0 -r Rp.nt.fasta_quant.sf -g Rp.v1_alt.fasta --dna Rp.nt.fasta -o test


Requires:
samtools 1.2 or later for Bam file
Biopython
NUmpy


Options:
  -h, --help            show this help message and exit
  --gff=FILE            gff file for predicted genes.
  --LTG=FILE            LTG outfile. This is the output generated  from the
                        Lateral_gene_transfer_prediction_tool
  --dna=FILE            predicted cds nucleotide genes for AT content stats
  -g FILE, --genome=FILE
                        genome.fasta
  -s SD_NUMBERS         the number of stadard deviations away from the mean
                        for identifying genes  that differ from normal AT
                        content. default=0
  -r RNASEQ, --rnaseq=RNASEQ
                        RNAseq expression profile for genes.  in format # Name
                        Length  TPM     NumReads  standard Sailfish output.
  -b BAM_FILE, --bam=BAM_FILE
                        bam file (sorted, indexed for samtools)  with genomic
                        reads mapped to geneome  this is used to see if HGT
                        genes have a different  genomic coverage compared to
                        other gene. Requires  samtools 1.2 or later
  -o OUT_FILE, --out_file=OUT_FILE
                        outfile to list the bad contigs


1) GENRATE bam file with genomic reads mapped to it:
How ever you want to do it, but sort and index your bam file
transrate --assembly genome.fasta --left genomic_reads.r1.fq.gz --right genomic_reads.r2.fq.gz --threads 12

BAM file is not need and can be run without it.  = much faster!!

2) GFF3
You may have to tidy and sort your GFF to a GFF3. Use GenomeTools
.. convert augustus.gft to gff3
.. gt gtf_to_gff3 -o test.gff3 -tidy augustus.gtf
or
.. gt gff3 -sort -tidy augustus.gff > formatted.gff3

3) LTG_LGT_candifates_AI_30plus.out:
This is the ouput from the Lateral_gene_transfer prediction tool. Precurser to this script.

4) RNAseq_coverage:
Agin, however you want to generate it. e.g.
transrate --assembly gene.cds --left rnaseq_r1.fq.gz --right rnaseq_r2.fq.gz --threads 12

5) Genome seq -g

6) cds of genes:
If you dont have it can use:
gffread *gff -g genome.fasta -x nt.fa -y aa.fa

BAD SCAFFOLDS??

The script will check to see if a contig is only made up of LTG/HGT predicted genes. If so, then this contig is suspect
and therefore should be considered as contimination.
Users are encouraged to used Blobplots of the genome assemblies before they get to this point.


