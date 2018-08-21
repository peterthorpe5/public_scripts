# README.py - `metapy_tools`
This repository is for development of ITS1-based diagnostic/profiling tools for the THAPBI Phyto-Threats project, funded by BBSRC.

In this directory, there are additional scripts for use with metapy

# USAGE:

SEE ``script.py -h`` on how to use the individual script.

``metapy.sh``   This is a shell script to run metapy on all the fastq or fastq.gz pairs in a directory. It doesnt care, compressed or not. It will work it out. 
You MUST fill in the variables virtual_machine, Working_directory and Raw_illumina_data. 
By deafult Swarm, Vsearch and Bowtie are run on the data. To run the other, either run metapy.py manually, or modify the ``generate_metapy.py`` script. The shell script will
output all the swarm and bowtie results in a folder obviously named. From the species list, it will identify where these were found in the samples. In the folders there is also,
total_number_of_assembled_sequences.txt and swarm/bowtie_percentage_cluster_with_db.txt. These numbers are not correct for vsearch, as the dereplication function
is done within vsearch not metapy, therefore these reads are not tracked. 
The total number of starting reads for each sample is in the sample results file. Grep, you can put these into one file, if you wish. 
Novel program clusters: This script, by default will go through all clusters ouputted by swarm, that did ot cluster with the databse, and BLAST search a representative of this cluster 
against nt database. This does NOT identify new species, but does identify other oomycetes and loose matches to Phytophora species. It does this using:
``bin/Novel_top_hits.py``  --min_cluster_size 65 --threads 16. So, only if a cluster is bigger than 65 it will BLAST this. Otherwise, the whole process will take to long.
Any smaller cluster is most likely not biologically relevant. Why this script? So, you dont have to call metapy,py 96 times. This will do it for you. The extra commands 
avavilble for metapy are not easily used from this script, so you will ave to modify generate_metapy.py for that. 


``generate_metapy.py``   This is a scrip to run metapy on all the reads in a folder. This is called from a shell script. 
This currently runs metapy in default mode. If you want to alter this, then you will have to manually add the commands into the
script. You may need to alter the path where it looks for metapy.py, line 33. Or put metapy.py in your path


``results_summary.py``  This is a scrip to summarise what has been found from all the swarm or bowtie results. Basically run from within the bowtie or swarm results 
directory. You must apply a threshold. e.g. -t 50, will only report results with more than 50 reads clustering with it. -o is the ouput name. If you run this within 
the novel_blast_results (discussed below), you need to add the --blast yes command. -m mismatches: only report results with -m or less mismatches from its database hit 
default is 5. 

``bin/BLAST_xml_parser.py``   The represnetative of the novel cluster is BLAST searched against nt (or a different db, you can give it that as an option, you may need to alter 
the path for nt anyaway), the ouput is an xml file. This script is called by bin/Novel_top_hits.py to filter the xml output file. 

``bin/Metapy_sentitivity_test.py``  This works out the specificity and snesitively based on our 4 groups of control mixes. 


shell_scripts:
==============
In here there are old shell scripts used to analyse various plates of data. These are now old but kept for reference. 

``Interogate_controls_all_folders.sh``   IN principle this script runs through all the control sequences and outputs the number of mismatches the sequences has against 
its correcsponding control sequence. Then the number of sequences with exactly 1,2,3,4 etc ... mismatches from its corresponding control sequences. Why? To see how 
much variation occurs after a single input sequence is put through the metabarcoding process. The script also outputs the coordinate of the mismatches so the user 
can see if specific regions are more prone to error. 

plot scripts:
=============
``identify_common_seq_plot_frequency.py``  For example te plot_files.sh uses this script to plot the frequency of unique sequences per sample. But late rin the script  
species are grouped together, for example all T30-4 samples, and plotted. if you add the -s yes command to this , the script will print to sdout the common seq for 
all input files (commor seperated). Therefore, all the T30-4 samples with the -s yes option will have 113 common ITS1 sequences. Are these real, or common errors. There is 
alot of evidence of multiple ITS1 regions within the genome.  "count_seq_abundance_from_related_samples.py" is an earlier version of this script. 

``Jellyfish_fq_to_table.py``  and  ``compare_jellyfish_mann_whitney.py``. _to_table.py script will count the orrurance of kmers are a user defined kmer lenght. This can be used 
to compare the sequences between two related samples. If they are the same then the kmer distribution will be the same. For example, single round PCR, versus nested. 
The nested is likely to introduce more errors and prone to more contamination. The compare_jellyfish_mann_whitney.py script will compare the kmer histograms and plot the histograms of the 
kmer frequencies. This showed that the single round versus nested PCR were statistically different from each other. 

