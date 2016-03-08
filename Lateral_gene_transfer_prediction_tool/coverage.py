
# module to get mapping coverage stats
# input is bam file

#imports
import numpy
import subprocess
import tempfile
from collections import deque
import datetime
import os


def get_total_coverage(bam_file, outfile):
    ## Run samtools idxstats (this get the coverage for all scaffolds:
    # assigne the outfile with the temp folder to keep thing more tidy
    oufile_dir = "./temp/"+outfile
    cmd = 'samtools idxstats "%s" > "%s"' % (bam_file, oufile_dir)
    #data was saved in idxstats_filename
    return_code = os.system(cmd)
    if return_code:
        clean_up()
        sys_exit("Return code %i from command:\n%s" % (return_code, cmd))
    # creat a dictioanry to hold all the total expression values for the scaffolds.
    overall_expression_dic = dict()
    with open(oufile_dir, "r") as handle:
        for line in handle:
            data = line.rstrip("\n").split("\t")
            scaffold = data[0]
            overall_expression_dic[scaffold] = [int(x) for x in data[1:]]
    #print ("scaffold Rp1: len,coverage,somethingelse = ", overall_expression_dic["Rp1"])
    #returns a dictionary: key[scaffold], vals = ['577', '274', '0'] len, reads_mapped, last_coloumn
    return overall_expression_dic

def average_standard_dev(positions):
    "function to return the avaerage for a list of number"
    the_mean = sum(positions) / float(len(positions))
    standard_dev = numpy.std(positions)
    return the_mean, standard_dev


def mean_coverage(coverage_array, slice_start, slice_end):
    selected_coverage = coverage_array[slice_start:slice_end]
    return mean(selected_coverage)


def stats(list_of_values):
    the_mean = sum(list_of_values) / float(len(list_of_values))
    #calc SD for AT for all genes
    standard_dev = numpy.std(list_of_values)
    return the_mean, standard_dev


#main function
def get_scaffold_coverage(genome, scaffold_to_gene_dict, gene_start_stop_dict, \
                          bam_file, overall_expression_dic,\
                          HGT_predicted_gene_set):
            #test one: is the a cds predcited for this, if not move on
    #overall_expression_dic = get_total_coverage(bam_file, "overall_scaffold_expression.txt")
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    list_of_coverage = []
    scaffold_mean_SD_cov_dict = dict()
    HGT_gene_to_genic_cov_dic = dict()
    for genome_record in SeqIO.parse(genome, "fasta"):
        # call samtools to get the depth per posititon for the seq of interest
        depth_filename = "./temp/depth.tmp"
        cmd = 'samtools depth -r "%s" "%s" > "%s"' % (genome_record.id, bam_file, depth_filename)
        #print("Running %s" % cmd)
        return_code = os.system(cmd)
        assert return_code == 0, "samtools says NO!! - something went wrong. Is your BAM file correct?"
        
        # assign zeros to all positions of the scaffold, as samtool does no report zeros
        all_coverage = [0]*len(genome_record.seq)
        for line in open(depth_filename):
            ref, possition, coverage = line.rstrip("\n").split("\t")
            list_of_coverage.append(coverage)
            possition = int(possition) - 1
            
            #assign the correct coverage to the relavent postiton, thus removing the zero value.
            #Or if there is no info, it remains as a zero.
            all_coverage[possition] = int(coverage)
            #Now have all the coverage information for this scaffold in all_coverage
            #get a list of gene on the current scaffold
            if scaffold_to_gene_dict[genome_record.id]:
                genes = scaffold_to_gene_dict[genome_record.id]
            else:
                continue              
                
            for i in genes:
                #is any of those gene in the HGT predicted list
                if i in HGT_predicted_gene_set:
                    # get its coordinates
                    start, stop = gene_start_stop_dict[i]
                    # get the average coverage of this region.
                    coverage_of_interest = mean_coverage(all_coverage,start,stop)
                    HGT_gene_to_genic_cov_dic[i]= coverage_of_interest
                        
            the_mean, standard_dev = average_standard_dev(all_coverage)
            data_formatted = "%f\t%f" %(the_mean, standard_dev)
            scaffold_mean_SD_cov_dict[genome_record.id] = data_formatted

    mean_genomic_cov, standard_dev_genomic_cov = average_standard_dev(list_of_coverage)

    return HGT_gene_to_genic_cov_dic, scaffold_mean_SD_cov_dict, mean_genomic_cov, standard_dev_genomic_cov


    























        
