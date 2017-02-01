READme info ITS_region_genomic_coverage
======================================================
author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

basic usage:

``./P.X.sh``

some variables in here need to be set by the user, if this piplien isused for other purposes.
These shells call the following script. (All regions in this file between ####
need to be filled in by the user).

``./ITS_genomic_coverage_ratio_finding.sh`` 

why?: We are interested in the raio of coverage of all identified
ITS regions in a genome versus "normal" genes - whatever they are.
Normal genes, as these cant be defined, will be treated as all other
genes, except RNA annotated genes (if annotated in GFF).

How?: Illumina reads are mapped back to the genome. Using the GFF3
gene coordingates, the coverage in obtained. A GFF3 file is generated
based on the BLAST results (ITS versus genome). The coordinates 
are used to obtain the genomic read coverage for these ITS regions.

Results?: The mean number of reads the map to "normal genes" versus 
ITS can be compared. This ratio give us an indication of the total
number of ITS regions that should be there, assuming genomic read 
coverage is normal and representative of gene number.


``Requires:``
	1) Internet access for $ wget
	2) Trimmomatic (read quality trimming)
	3) BLAST
	4) python
	5) bowtie2
	6) samtools
	7) bedtools   *
	8) BUSCO version 1.1b - currently does not run with version 1.2
		
*(for poorly formamted GFF files - try reformmatting with GenomeTools)

BUSCO is used to predict the "EOG" genes from the genomes. WHY? 
basically very few of them have publically availble gene models.
To get a backgroud gene count, and one that is supposed to be, as 
far as possible, based on single copy genes - BUSCO is used. 

You will have to prepare the data for BUSCO before hand by dowloading
http://busco.ezlab.org/files/eukaryota_buscos.tar.gz
.. put the path to this, once decompressed in the shell scripts as
required.


