#!/usr/bin/env python
#title: Parse clusters and find the phy spcies
#in the cluster
#author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

#imports
import os
import sys
from optparse import OptionParser
import datetime
import os
from sys import stdin,argv
import subprocess

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator



#####################################################################


def stop_err( msg ):
    """stop function"""
    sys.stderr.write("%s\n" % msg)
    sys.exit(1)

try:
    # New in Python 3.4
    from statistics import mean
except ImportError:
    def mean(list_of_values):
        """Calculate the mean average of a list of numbers."""
        # Quick and dirty, assumes already a list not an interator
        # so don't have to worry about getting the divisor.
        # Explicit float(...) to allow for Python 2 division.
        return sum(list_of_values) / float(len(list_of_values))

assert mean([1,2,3,4,5]) == 3

def count_reads(in_fastq):
    """function count the number of reads"""
    #open the fastq file
    in_file = open(in_fastq)
    # iterate through the fastq file
    total_reads = 0
    for (title, sequence, quality) in FastqGeneralIterator(in_file):
        total_reads = total_reads+1
    in_file.close()
    return total_reads

def get_fasta_stats(fasta):
    """function to get stats on a given fasta file"""
    with open(fasta, 'r') as seq:
        sizes = [len(record) for record in SeqIO.parse(seq, 'fasta')]
    min_contig = min(sizes)
    max_contig = max(sizes)
    avg_contig = mean(sizes)
    num_contig = len(sizes)
    return min_contig, max_contig, avg_contig, num_contig

def parse_line(line):
    """finction to parse line"""
    if not line.strip():
        return False #if the last line is blank
    if line.startswith("#"):
        return False
    #if not "Phytophthora" or "P." in line:
        #continue
    line = line.rstrip()
    out_put_str = ""
    cluster_line = line.rstrip("\n").split("\t")
    # are the database Phy singletons? - then we wont be interested
    number_of_reads_hitting_species = len(cluster_line)
    #set up some blank variables
    species = ""

    #database phy counter
    phy_count = 0
    return cluster_line, number_of_reads_hitting_species, species, phy_count

def make_blastdb(infile):
    """function to make a nucl blastdb"""
    cmd = 'makeblastdb -in %s -dbtype nucl' % (infile)
    #print("Running %s" % cmd)
    return_code = os.system(cmd)
    assert return_code == 0, "blast makedb says NO!! - something went wrong. Is your fa file correct?"
    
def run_blastn(infile, cluster_count):
    """function to run blastn using"""
    cmd = 'blastn -db %s -query %s -outfmt 6 -out %s.tab' % (infile, infile, infile.split(".fasta")[0])
    #print("Running %s" % cmd)
    return_code = os.system(cmd)
    assert return_code == 0, "blastn says NO!! - something went wrong. Is your fa file correct?"
    cmds = 'rm %s.n*' % (infile)
    #print("Running %s" % cmds)
    return_code = os.system(cmds)

def align_cluster(infile):
    """function to muscle align a cluster"""
    cmd = 'muscle -in %s -out %s_TEMP' % (infile, infile)
    print("Running %s" % cmd)
    return_code = os.system(cmd)
    assert return_code == 0, "muscle says NO!! - is muscle in your PATH?"
    #wait command forces Linux to wait until all command are complete
    cmd_wait = 'wait'
    print("Running %s" % cmd_wait)
    return_code = os.system(cmd_wait)
    #refine take longer but is the most accurate
    cmd2 = 'muscle -in %s_TEMP -out %s_aligned.fasta -refine' % (infile, infile.split(".fasta")[0])
    print("Running %s" % cmd2)
    return_code = os.system(cmd2)
    assert return_code == 0, "muscle -refine says NO!! - is muscle in your PATH, does this version allow this parameter?"
    #wait command forces Linux to wait until all command are complete
    cmd_wait = 'wait'
    print("Running %s" % cmd_wait)
    return_code = os.system(cmd_wait)

    cmd3 = 'rm %s*_TEMP' % (infile)
    print("Running %s" % cmd3)
    return_code = os.system(cmd3)
    assert return_code == 0, "something went wrond deleting the temp alignment file!" 

def parse_tab_file_get_clusters(fasta, all_fasta, in_file,min_novel_cluster_threshold, \
                            read_prefix, show_me_the_reads, right_total_reads, working_directory,\
                            v, blast, align, out_file):
    """#script to open up a tab separeted clustering output and identify the
    species in the clustering"""
    #print coded_name_to_species_dict

    #index the all_seq fasta.
    all_sequences =  SeqIO.index(all_fasta, "fasta")
    #call the function to get fasta stats
    min_contig, max_contig, avg_contig, num_contig = get_fasta_stats(fasta)
    read_prefix = str(read_prefix)

    cluster_file = open (in_file, "r")
    summary_out_file = open(out_file, "w")

    ITS_hitting_phy_count = 0
    #title for the results file
    title = "#cluster_number\tspecies\tnumber_of_reads_hitting_species\treads_that_hit_species\n"
    summary_out_file.write(title)
    cluster_number = 0
    for line in cluster_file:
        cluster_number +=1
        #call function to parse line
        if not parse_line(line):
            continue
        reads_of_interest = ""
        cluster_line, number_of_reads_hitting_species, \
                      species, phy_count = parse_line(line)
        # basically not a singleton cluster
        if len(cluster_line) >1:        
            #open a file to put seq into
            out_name = "%s/clusters_d%d/cluster%d_len%d.fasta" %(working_directory, v, cluster_number, len(cluster_line))
            cluster_fasta = open(out_name, "w")

        for member in cluster_line:
            member = member.rstrip()
            try:
                seq_record = all_sequences[member]
                seq_record.description = ""           
            except:
                ValueError
                stop_err ("""missing sequence %s. This should not be seen. Have you passed me the correct
                        fasta file?\n\n""" % member)# break the program here
            try:
                SeqIO.write(seq_record, cluster_fasta, "fasta")
            except:
                ValueError
                #print ("singlton cluster - no file generated for cluster %d" %cluster_number)
            #is a memeber of the database in this cluster?
            if "Phytophthora" in member or "P." in member or "VHS" in member:
                # yes we are interested in this cluster
                phy_count = phy_count+1
                #print member
                species = species+member+" "
                #remove all the memebr which are database memebers
                number_of_reads_hitting_species = number_of_reads_hitting_species-1
        if number_of_reads_hitting_species >=1 and phy_count >0:
            #after the line, are there more memeber that are not database members?
            for member in cluster_line:
                if "Phytophthora" in member or "P." in member or "VHS" in member:
                    continue
                if not member.startswith(read_prefix):
                    continue
                # we do not add the reads that hit the 
                reads_of_interest = reads_of_interest+member+" "
  
            ITS_hitting_phy_count = ITS_hitting_phy_count+number_of_reads_hitting_species
            #format the data
            if show_me_the_reads:
                data_output = "%d\t%s\t%d\t%s\n" %(cluster_number, species.rstrip(), \
                                                   number_of_reads_hitting_species,\
                                                   reads_of_interest)
            else:
                data_output = "%d\t%s\t%d\n" %(cluster_number, species.rstrip(), \
                                               number_of_reads_hitting_species)
                
            #write out the data
            summary_out_file.write(data_output)
            phyto_read = (float(ITS_hitting_phy_count)/ num_contig)*100
        else:
            if len(cluster_line) > min_novel_cluster_threshold:
                #print "cluster len", len(cluster_line), "threshold", min_novel_cluster_threshold, cluster_line
                #print ("novel cluster cluster %d" %(cluster_number))
                #open a file to put seq into
                out_novel = "%s/novel_d%d/novel%d_len%d.fasta" %(working_directory, v, cluster_number, len(cluster_line))
                novel_fasta = open(out_novel, "w")
                for member in cluster_line:
                    member = member.rstrip()
                    seq_record = all_sequences[member]
                    seq_record.description = ""           
                    SeqIO.write(seq_record, novel_fasta, "fasta")
                    
        #close the open files, makeblastdb run blastn on individual clusters to get % identify
        try:
            cluster_fasta.close()
            #call functions to run blast on the cluster
            if blast:
                make_blastdb(out_name)
                #run blast and remove db files
                run_blastn(out_name, cluster_number)
            if align:
                #align the cluster
                #print ("I am about to align the cluster: %s" % out_name)
                align_cluster(out_name)
 
        except:
            ValueError #singleton cluster - no file generated
        try:
            novel_fasta.close()
            #call functions to run blast on the cluster
            if blast:
                make_blastdb(out_novel)
                #run blast and remove db files
                run_blastn(out_novel, cluster_number)
            if align:
                #align the cluster
                #print ("I am about to align the cluster: %s" % out_novel)
                align_cluster(out_novel)
            
        except:
            ValueError #no novel files to close

    fasta_file_summary = """#Fasta file assembly summary: #min_contig = %d max_contig = %d avg_contig = %d
#Total number of assemblerd sequences = %d
#number of reads clustering with Phy = %d
#number of starting reads = %d
#percent of reads clustering with Phyto = %.2f""" %(min_contig, \
                                            max_contig, avg_contig, num_contig, ITS_hitting_phy_count,\
                                            right_total_reads, phyto_read)
    summary_out_file.write(fasta_file_summary)
    cluster_file.close()
    summary_out_file.close()
    return True

#############################################################################

#to run the script       

usage = """usage :

this scripts to summarise what clusters with what Phy species. The script also
dumps one file per cluster, BLASTN searches these back against itself
and aligns it.

WARNING: THIS WILL PRODUCE 10 000S OF THOUSANDS OF FILES WHICH WILL SLOW YOU COMPUTER
DOWN, AND MAKE YOU UNPOPULAR WITH SYS ADMIN. YOU ASKED FOR THIS FEATURE!!! The BLAST and alignment
also makes it very slow. 


$ python parse_clusters.py -i clustering_outfile_already_decoded_from_temp_names -n 3 -a all_seqeucnes.fasta
--read_prefix M01157 -o summarise_clusters.out

command line option


"""

parser = OptionParser(usage=usage)

parser.add_option("-f","--fasta", dest="fasta", default=None,
                  help="assembled fasta file to get stats on",
                  metavar="FILE")

parser.add_option("-a","--all_fasta", dest="all_fasta", default=None,
                  help="both the assembled fasta and databse fasta."
                  " this is to get the sequences from from for the clusters folder",
                  metavar="FILE")

parser.add_option("-n","--min_novel_cluster_threshold", dest="min_novel_cluster_threshold", default=4,
                  help="min_novel_cluster_threshold to determine is this is a"
                  " novel clusters to output. Default = 4 ")

parser.add_option("-l", "--left", dest="left", default="temp_not_trimmedr1.fq",
                  help="left reads unzipped fq file. needed to get count of reads. ",
                  metavar="FILE")

parser.add_option("--Name_of_project", dest="Name_of_project", default=None,
                  help="name of project to make folders. ")

parser.add_option("-i","--in", dest="in_file", default=None,
                  help="clustering out file")

parser.add_option("-v","--difference", dest="v", default=None,
                  help="current swarm v value")

parser.add_option("-r", "--right", dest="right", default="temp_not_trimmedr2.fq",
                  help="right reads unzipped fq file. needed to get count of reads. ",
                  metavar="FILE")

parser.add_option("-s", "--show_me_the_reads", dest="show_me_the_reads", default=False,
                  help="show_me_the_reads in the output file for those that hit the Phy species."
                  " by default this is off, as the file could get very large. -s True if you want ...",
                  metavar="FILE")

parser.add_option("-t", "--threads", dest="threads", default="4",
                  help="number of threads for blast - threads ",
                  metavar="FILE")

parser.add_option("--read_prefix", dest="read_prefix", default=None,
                  help="read_prefix from the fq file. Only needs "
                  " the first few letter e.g. read_prefix  M01157 ",
                  metavar="FILE")

parser.add_option("--blast", dest="blast", default=False,
                  help="this option performs BLAST on itself. Turn off for faster results."
                  " To turn of --blast False ")

parser.add_option("--align", dest="align", default=False,
                  help="this option performs alignment on a cluster. Turn off for faster results."
                  " To turn off --align False")

parser.add_option("-o", "--out_prefix", dest="out_file", default="summarise_clusters.out",
                  help="prefix to the output filenames")


(options, args) = parser.parse_args()
#-f
fasta = options.fasta
# --all_fasta
all_fasta = options.all_fasta
#-i
in_file = options.in_file
# -l
left = options.left
# -r
right = options.right
# -s
show_me_the_reads = options.show_me_the_reads
# -t
threads = options.threads
# -n
min_novel_cluster_threshold = options.min_novel_cluster_threshold
# -o
out_file = options.out_file
# --blast
blast = options.blast
# --align
align = options.align
# -v
v = options.v
# --Name_of_project
Name_of_project = options.Name_of_project
#read_prefix
read_prefix = options.read_prefix



###################################################################
#start of program

#call the function to count the trimmed reads 
right_total_reads = count_reads(right)
left_total_reads = count_reads(left)

#sanity test.
assert right_total_reads == left_total_reads, """ \n\nthe total number of reads in the
left and right file do not match. Something has gone wrong before here! These should
be the same. Are you giving me the the correct files?\n\n"""

#make sure this is an int
min_novel_cluster_threshold = int(min_novel_cluster_threshold)
v = int(v)

###################################################################
#make folders to put the cluster fasta files. 
file_name = 'test.txt'
make_folder_list = ["clusters", "novel"]
working_dir = os.getcwd()
for name in make_folder_list:
    name = name+"_d%d" % (v)
    working_directory = os.path.join(working_dir, Name_of_project+"_results")    
    dest_dir = os.path.join(working_dir, Name_of_project+"_results", name)
    try:
        os.makedirs(dest_dir)
    except OSError:
        print ("folder already exists, I will write over what is in there!!")

###################################################################
#run the program
parse_tab_file_get_clusters(fasta, all_fasta, in_file,min_novel_cluster_threshold, \
                            read_prefix, show_me_the_reads, right_total_reads, \
                            working_directory, v, blast, align, out_file)

print "done"
