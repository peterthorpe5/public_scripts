#!/usr/bin/env python

########################################################################################################
#####: TITLE: script to get the upstream regions of genes of interest ##################################
#script will return upt to the gene if the full length falls within that gene.
#also, script will return reverse complemnet of negative strand coded genes.

#author: Peter Thorpe September 2015. The James Hutton Insitute, Dundee, UK.

"""
script to return a user defined threshoold upstream number of nucleotides of genes of interest
from the genome seq using a modified GFF output. 
"""

#biopython imports
from Bio.Seq import Seq
from Bio import SeqIO
import time
#os imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser


########################################################################################################
def wanted_genes(genes_file):
    "function to return a list of wanted genes from file"
    wanted = open(genes_file, "r")
    names = wanted.readlines()
    wanted_data = [line.rstrip() for line in names
              if line.strip() != "" if not line.startswith("#")]
    wanted.close()
    wanted_set = set([])
    for i in wanted_data:
        i = i.replace("id=", "")
        i = i.replace(".t1", "")
        wanted_set.add(i.split(";")[0])
    return wanted_set

def index_gene_scaffold_coordinates(coordinate_file):
    "function to return dictionary genes and coordinates without directions"
    coordinate_dict = dict()
    data = open(coordinate_file, "r")
    genes_coordinate = data.readlines()
    genes_coordinate = [line.replace("ID=", "").rstrip() for line in genes_coordinate
              if line.rstrip() != "" if not line.startswith("#")]
    gene_list = []
    data.close()
    for gff_info in genes_coordinate:
        #print gff_info.split("\t")[4]
        assert len(line.split("\t")) ==5 ,"GFF/ altered file fields wrong length should be 5"
        gene=gff_info.split("\t")[4]
        if ";" in gene:
            gene = gene.split(";")[0]
        #check each gene only represented once
        assert gene not in gene_list, "duplicate genes found. Reformat file -C file. Problem gene is: %s" % (gene)
        gene_list.append(gene)
        if gene in coordinate_dict.values():
            print "repeated line in gff sub file"
            continue
        else:
            scaffold_cordinates = gff_info.split("\t")[:]
            coordinate_dict[gene] = scaffold_cordinates
    #print coordinate_dict
    return coordinate_dict

  
def iterate_through_coordinate_dictionary(coordinate_dict,gene_name, scaffold, coordinates):
    "check to see if the scaffold and new coordinate hits a predicted gene"
    for gene, vals in coordinate_dict.items():
        #find the genes on the same scaffold
        if scaffold in vals[0]:
            dictionary_scaffold = vals[0]
            gene = vals[4]
            #if its is the same gene as the stop
            if gene_name ==gene:
                continue
            if scaffold != dictionary_scaffold:
                continue
                
            else:
                #debugging comment due to Roman numeral scaffold name being "within" eachother
                #print "scaffold = ", scaffold, "acutally looking at", vals[0]

                #cureent gene start
                start = int(vals[1])
                #cureent gene stop
                stop = int(vals[2])
                coordinates = int(coordinates)
                direction = vals[3]
                #print "coordinates=%s, start=%s, stop=%s , gene_query=%s, gene_GFF=%s" %(coordinates, start, stop, gene_name, gene)
                #basically does the coordinate fall in the current coordinate for a gene
                if coordinates >= start and coordinates <= stop:
                    print "\n%s gene upstream coordinate falls in the genic regions of %s on %s scaffold\n" % (gene_name,gene,scaffold)
                    data = coordinate_dict[gene_name]
                    if "+" in data:
                        #+ coding gene, upstream will be returned as the end (stop) of the preceding gene
                        return stop
                    else:
                        #- coding gene, upstream will be returned as the begining (start) (stop) of the proceding gene
                        return start
    return False


    
def parse_through_gene_coordinates(coordinate_file):
    "function to retunr genes and coordinates with directions"
    data = open(coordinate_file, "r")
    genes_coordinate = data.readlines()
    genes_coordinate = [line.replace("ID=", "").rstrip() for line in genes_coordinate
              if line.strip() != "" if not line.startswith("#")]
    data.close()
    return genes_coordinate


def seq_getter(coordinate_file, genome_sequence, upstream, genes_file, outfile, user_defined_genic=0):
    """this is a function returns the upstream regions of the list of genes of interest
     - a user defined threshold for the number of nucleotides upstream is also used"""
    f= open(outfile, 'w')
    threshold = int(upstream)

    wanted = wanted_genes(genes_file)
    assert len(wanted) == len(set(wanted)), "duplicates in wanted list!"

    Genome_sequence =  SeqIO.index(genome_sequence, "fasta")
    Genome_sequence_time=time.time()
    #print 'import genome seq file took, %.3f' %(Genome_sequence_time - start_time)
    coordinates = parse_through_gene_coordinates(coordinate_file)
    #[gene] = scaffold, cordinates
    cordinate_dictionary = index_gene_scaffold_coordinates(coordinate_file)
    #we need to get the gene start and stop information
    #assign the user defined number of genic nucleotides to return to int
    user_defined_genic = int(user_defined_genic)
    for gene_name in wanted:
        coordinate_info = cordinate_dictionary[gene_name]
        #assert len(coordinate_full_info.split("\t"))==5, "coordinate infor not prperly formatted.! I want: gene	start	stop	+/-	ID=scaffold_number"
        contig = coordinate_info[0]
        location_start = coordinate_info[1]
        location_stop = coordinate_info[2]
        direction_of_coding = coordinate_info[3]
        #get the altered start (+) strand, and stop (-) strand.
        #if the user want genic regions, user_defined_genic will bring back X amount of the gene
        DNA_start =int(location_start)
        DNA_stop = int(location_stop)
        
        assert DNA_start<DNA_stop, "are you sure the start stop coordinates are correct?"
        upstream = (DNA_start)-threshold
        if upstream<2:
            print "\nWARNING: %s upstream region may fall off start of contigs %s. Check this\n" % (gene_name, contig)
            
    
        # yes DNA_stop - this is how it is coded in GFF files
        neagtive_strand_upstream = (DNA_stop)+threshold

        Genome_seq_record = Genome_sequence[contig]
        length_of_contig = len(Genome_seq_record.seq)
        DNA_region_of_interest = Genome_seq_record.seq

        if neagtive_strand_upstream > length_of_contig:
            print "WARNING: %s upstream region may fall off end of contig %s. Check this\n" % (gene_name,contig)
            
                   #test if the "upstream region hits a gene


        #slice up the contig for the region of interest [threshold upstream: genestart]or negative [geneend: genestart plus threshold
        DNA_region_of_interest_upstream_positive = DNA_region_of_interest[upstream:(DNA_start+(user_defined_genic-1))]

        
        #for negative strand - we reverse complement it
        DNA_region_of_interest_negative_upstream2 = DNA_region_of_interest[(DNA_stop-user_defined_genic):(neagtive_strand_upstream-1)]
        DNA_region_of_interest_negative_upstream = DNA_region_of_interest_negative_upstream2.reverse_complement()

        if "-" in direction_of_coding:
            new_start = iterate_through_coordinate_dictionary(cordinate_dictionary, gene_name, contig, neagtive_strand_upstream)
            if new_start:
                new_start= int(new_start)
                print "this %s query request is going to return a region of a gene" %(gene_name)
                print "I am only going to return up to %s  NNNEEWWW SSTTAARRTT" %(new_start)
                
                DNA_region_of_interest_negative_upstream2 = DNA_region_of_interest[(DNA_stop-(user_defined_genic)):new_start]
                DNA_region_of_interest_negative_upstream = DNA_region_of_interest_negative_upstream2.reverse_complement()
                if len(DNA_region_of_interest_negative_upstream)>50:

                    print >> f, '>%s\t|%s\t[%s:%s]%sbp_upstream - strand\n%s' % \
                            (gene_name,contig,(DNA_stop-user_defined_genic),
                             str(new_start), threshold,
                             DNA_region_of_interest_negative_upstream)
                    if "NNNNN" in DNA_region_of_interest_negative_upstream:
                        print "NNN found in %s" %(gene_name)
            else:
            #(-) ... downstream is (DNA_stop) + threshold"
                if len(DNA_region_of_interest_negative_upstream) >50:
                    print >> f, '>%s\t|%s\t[%s:%s]%sbp_upstream - strand\n%s' % \
                            (gene_name,contig,(DNA_stop-user_defined_genic),
                             neagtive_strand_upstream, threshold,
                             DNA_region_of_interest_negative_upstream)
                    if "NNNNN" in DNA_region_of_interest_negative_upstream:
                        print "NNN found in %s" %(gene_name)
            
        if "+" in direction_of_coding:
            #test if the "upstream region hits a gene
            new_region_to_return = iterate_through_coordinate_dictionary(cordinate_dictionary, gene_name,contig, upstream)
            #print "newregion ", new_region_to_return
            if new_region_to_return:
                new_region_to_return = int(new_region_to_return)
                #print "old start %s" %DNA_region_of_interest_upstream_positive
                #print "this %s query request is going to reurn a region of a gene" %(gene_name)
                DNA_region_of_interest_upstream_positive = DNA_region_of_interest[new_region_to_return:(DNA_start+(user_defined_genic-1))]
                #print "new start %s" %DNA_region_of_interest_upstream_positive
                if len(DNA_region_of_interest_upstream_positive) >50:
                    print >> f, '>%s\t|%s\t[%s:%s]%sbp_upstream + strand\n%s' % \
                        (gene_name, contig,\
                         str(new_region_to_return), (DNA_start+user_defined_genic), threshold,
                         DNA_region_of_interest_upstream_positive)
                    if "NNNNN" in DNA_region_of_interest_upstream_positive:
                        print "NNN found in %s" %(gene_name)

            else:
            #(+) ... upstream is (DNA_start) - threshold"
                if len(DNA_region_of_interest_upstream_positive)>50:
                    print >> f, '>%s\t|%s\t[%s:%s]%sbp_upstream + strand\n%s' % \
                        (gene_name, contig,\
                         upstream, (DNA_start+user_defined_genic), threshold,
                         DNA_region_of_interest_upstream_positive)
                    if "NNNNN" in DNA_region_of_interest_upstream_positive:
                        print "NNN found in %s" %(gene_name)
                          

##              #alterantive print methods if someone wnats properly formatted fasta files
##                #print >> f, seq_record.format("fasta")
##                #SeqIO.write call (with line wrapping), faster but still a bit slow
##            #SeqIO.write(seq_record, f, "fasta")
   
    f.close()
    return True



#####################################################################################################
start_time=time.time()
###################################################################################################


if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.3"
    sys.exit(0)


usage = """Use as follows:

$ python get_upstream_regions_using_gene_coordinates_GFF_format.py --coordinates coordinate_file.fasta -g genome_sequence.fasta -upstream <int> number of nucleotides upstream of strat of gene to return e.g.  -u 1000
-z user_defined_genic (how much of the gene to return) -o outfile_name

Requirements:
python 2.7 and biopyton.  If using python 3.x, the user must alter the print to file statements accordingly (print >> is 2.X only). 

This will return (--upstream number) of nucleotides to the start of your genes(s) of interest (-g) gene_file using data from (-c). Gene file can either be space, tab or \n separated..

The coordinate file can be generated using a GFF3 file and a linux command line using:

grep "gene" name.gff3 | grep -v "#" | cut -f1,4,5,7,9 > format_for_py_script.out


yeilding this reulting file:

scaffold	start	stop	strand(+/-)	ID=gene_number
GROS_00001	2195	3076	-	ID=GROS_g00002
GROS_00001	8583	10515	+	ID=GROS_g00005.....

The script will check that the start is always less than the end. GFF file should have starts < stop irrespective of the coding direction

To get all the genes in the file do:

cut -f5 format_for_py_script.out > all_gene_names.out

MORE help / example:

This is an example I ran for G. pallida:

python ~/misc_python/up_stream_genomic_regions/get_upstream_regions_using_gene_coordinates_GFF_format_plus_gene_region.py -c format_for_py_script.out -g Gpal.v1.0.fas -f all_gene_names.out -u 225 -z -125 -o Gp.all_gene_names.out_225up_125genic.fasta > warning_all_gene_names.out_225up_125genic.out


This got (-u 225) 225 bp upstream and (-z 125) 125bp into the current gene for all the genes in (-f) all_gene_names.out. By default -z is zero. So you dont need to specify this, unless you specifically want a piece of the current gene being searched for.
"""

parser = OptionParser(usage=usage)

parser.add_option("-c", "--coordinates",dest="coordinate_file", default="format_for_py_script.out",
                  help="NOTE: coordinate_file can generate using linux command line of "
                  "GFF file:  grep 'gene' name.gff3 | grep -v '#' | cut -f1,4,5,7,9 > format_for_py_script.out ."
                  "Default = format_for_py_script.out")

parser.add_option("-g", "--genome", dest="genome_sequence", default=None,
                  help="genome_sequence.fasta  -  this has to be the file used to generate the gene models/GFF file")

parser.add_option("-f", "--gene_names", dest="genes_file", default=None,
                  help="a file with a list of gene names to get the upstream regions for")

parser.add_option("-u", "--upstream", dest="upstream", default="1000",
                  help="the amount of nucleotide upstream of the gene start, taking into account gene directions, to return in the outfile"
                  "by default this will not return sequences of 50bp or less. If you require these alter"
                  "lines 204 and 214")

parser.add_option("-z", "--user_defined_genic", dest="user_defined_genic", default="0",
                  help="the number of nucleotides from within the gene to return, default is 0")

parser.add_option("-o", "--output", dest="out_file", default="upstream_of_genes.fasta",
                  help="Output filename (fasta file)",
                  metavar="FILE")


(options, args) = parser.parse_args()

coordinate_file = options.coordinate_file
genome_sequence = options.genome_sequence
upstream = options.upstream
genes_file = options.genes_file
outfile = options.out_file
user_defined_genic = options.user_defined_genic


#if len(args) < 1:
    #stop_err("Expects no argument, one input filename")

if not os.path.isfile(coordinate_file):
    sys_exit("Input coordinate_file file not found: %s" % coordinate_file)

if not os.path.isfile(genome_sequence):
    sys_exit("Input genome_sequence file not found: %s" % genome_sequence)

if not os.path.isfile(genes_file):
    sys_exit("Input genes_file file not found: %s" % genes_file)

seq_getter(coordinate_file, genome_sequence, upstream, genes_file, outfile, user_defined_genic)



end_time=time.time()
print 'that took, %.3f' %(end_time - start_time)

#print 'done'
