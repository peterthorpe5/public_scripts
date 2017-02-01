READme info Phyt_ITS_identifying_pipeline
======================================================
IMPORTANT:
THIS IS STILL IN TESTING PHASE. THIS IS A DRAFT PIPELINE AND WILL BE CONVERTED
TO PYTHON ONCE TESTING IS COMPLETE.

author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

basic usage:

run an individual shell script, which acts as a config file, that then call the main 
script below 
# config file - see bottom of here for how to generate this file. 
``./Phy_CONFIG_file.sh``

#which calls
``./Phyt_ITS_identify_pipline.sh``

or

``./Phyt_ITS_identify_pipline_with_error_correction.sh``

Important parameters to be set:
some variables in here need to be set by the user.
see top section marked: ## THESE VARIABLE NEED TO BE FILLED IN BY USER !!!
until: # NOTHING TO BE FILLED IN BY USER FROM HERE!



why?: We are interested in identifying Phyt species based on the ITS
regions. Phy species could become more previlent in the UK due to climate change and
the importation of plants.

How?: Following plant sampling, DNA extraction, primers amplify the ITS region.
PCR amplified products are Illumina sequenced. The reads are QC, trmimmed, QC,
assembled together (reads are overlapping) and clustered with a "known" database of Phyt
ITS regions. The resulting clusters are used to identify the Phy species present.

Results?: NOT FINISHED. Basically, the clustering is one cluster per line, tab separated.



``Requires:``
	1) Trimmomatic (read quality trimming) http://www.usadellab.org/cms/?page=trimmomatic
	2) Fastqc http://www.bioinformatics.babraham.ac.uk/projects/download.html
	3) python2.7 or greater. Python 3.5 is recommended. https://www.python.org/downloads/
	Biopython is also required:  http://biopython.org/wiki/Download
	4) PEAR * (assemble overlapping reads) https://github.com/xflouris/PEAR  http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.10-bin-64.tar.gz 
	5) Swarm (clustering) https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf
	6) To draw graphs: No essential but Numpy, Matplotlib, Scipy... pip install module_name
	7) BLAST. This is used to get the % identity of each cluster Please add to your PATH https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
	8) muscle. Alignmnet tool to align the clusters. Please rename and add to your PATH as muscle.  http://www.drive5.com/muscle/
	9) 	BLESS error correction (difficult to install) - Not essential. NOT CURRENTLY USED
	10) SPAdes alternative error correction. Much easier to install . http://bioinf.spbau.ru/en/content/spades-download-0 

* for PEAR:	
https://github.com/xflouris/PEAR/blob/master/bin/pear-0.9.5-bin-64  - precompiled and includes [+bzlib +zlib] included which are a pain to install.

https://github.com/xflouris/PEAR/tree/master/bin

generate config file
=====================
Either use the example config file as a template and fill in your experiment details.
Dont forget to make the file executable:
$chmod a+x names_of_config_file.sh``


or

use the python script:

``python generating_config_files.py``


Why these tools?
================
Pipeline: A species identification pipeline given Illumina data
A number of pipelines exist that have been developed to identify species based of OTUs, however mostly for fungi: FHiTINGS (Dannemiller, et al. 2014), Mycofier (Delgado-Serrano, et al. 2016), PIPITS (Gweon, et al. 2015). However, to my knowledge only one pipeline exists for use with Phytophthora (Scibetta, et al. 2012).
Therefore a pipeline needs to be developed for Phytophthora species, which utilises a Phytophthora specific ITS1 database. 
How
A number of tools exist for each step. The following listed are the current “best” tools for the job, in my opinion. Figure 1 represents the major steps required for the identification of species when given Illumina data. A known ITS1 database is required (OTUs). 


There are a number of different tools which can be used at each step. Arguably assembly and clustering tools can have the biggest impact. For assembly PEAR (Zhang, et al. 2014) was chosen. PEAR is very fast and applies statistical tests to minimise false positive results. For small overlaps it was reported that FLASH (Magoč and Salzberg 2011) often fails to merge such reads (Zhang, et al. 2014), therefore PEAR was chosen. 
For clustering, a number of tools exist. The data analysed in (Scibetta, et al. 2012) was clustered with BlastClust and OTUs picked with Qimme. Other clustering tools exist such as CD-HIT-454, DNAclust, Swarm and Usearch. Tools such as CD-HIT-454 and Usearch have inherent problem: They are biased to the order the dataset is presented in. Thus results can alter depending on the order of the dataset. Swarm (Mahé, et al. 2015) is an amplicon clustering tool that has been developed which aims to overcome the limitations of other clustering tools. To summarise the clusters, custom Python scripts are made. 


SWARM INFO:

-d, --differences "zero or positive integer"
maximum number of differences allowed between two amplicons, meaning
that two amplicons will be grouped if they have "Iinteger" (or
less) differences. This is swarm's most important parameter. The
number of differences is calculated as the number of mismatches
(substitutions, insertions or deletions) between the two amplicons
once the optimal pairwise global alignment has been found (see
"pairwise alignment advanced options" to influence that step). Any
Iinteger between 0 and 256 can be used, but high values
will decrease the taxonomical resolution of swarm
results. Commonly used values are 1, 2 or 3, rarely
higher. When using  0, swarm will output results
corresponding to a strict dereplication of the dataset, i.e. merging
identical amplicons. Warning, whatever the value, swarm
requires fasta entries to present abundance values. Default number of
differences is 1.



