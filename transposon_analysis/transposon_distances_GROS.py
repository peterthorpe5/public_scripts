#title: scrip to compare the distances of predicted transposons
# to all genes.

#why? transposon were found to be close to effector islands and HGT
# gene by using these scripts.

#author: Peter Thorpe September 2015. The James Hutton Insitute, Dundee, UK.

#imports
import os
import sys
from optparse import OptionParser
import datetime

####################################################################################
def split_gff_gene_names(gene_info):
    "function to hopefully return just the gene name"
    #print gene_info
    try:
        gene = gene_info.split("ID=")[1]
    except:
        gene = gene_info
    try:
        gene = gene_info.split("Parent=")[1]
    except:
        gene = gene
    gene = gene.split(".gene")[0]
    gene = gene.split(".exon")[0]
    gene = gene.split(";")[0]
    gene = gene.split(";")[0]
    gene = gene.split(".CDS")[0]
    #some data set require this to be uncommented
    #gene = gene.split(".t")[0]
    gene = gene.rstrip("\n")
    return gene

def parse_gff(gff):
    "function to parse GFF and produce a scaffold_to_gene_dict"
    f_in = open(gff, "r")
    # data to output of function
    scaffold_to_gene_dict = dict()
    gene_to_exon_count = dict()
    gene_start_stop_dict = dict()
    count = 1
    gene_name = ""
    #iterate through gff
    for line in f_in:
        if line.startswith("#"):
            continue
        scaffold,aug,cds_type,start,stop,e,f,g,gene_info = line.split("\t")
        gene = split_gff_gene_names(gene_info)
        #print gene
        gene_names = gene_names+"\n"+gene

        #scaffold_to_gene_dict
        if line.split("\t")[2] == "gene":
            if not scaffold in scaffold_to_gene_dict:
                scaffold_to_gene_dict[scaffold]=[gene]
            else:
                scaffold_to_gene_dict[scaffold].append(gene)
            start_stop_formatted = "%s\t%d\t%d" %(scaffold, int(start), int(stop))
            gene_start_stop_dict[gene] = start_stop_formatted

        #gene_to_exon_count
        if line.split("\t")[2] == "exon":
            if not gene in gene_to_exon_count:
                gene_to_exon_count[gene] = count
            else:
                gene_to_exon_count[gene]= count +1
    #print scaffold_to_gene_dict
    f_in.close()
    return gene_names, scaffold_to_gene_dict, gene_start_stop_dict

def LTG_file(LTG_file):
    """function to parse LTG prediction
    get a set of names, and a gene_to_comment_dict.
    The majority of this funcion was used for something else."""
    # in file is the output of the LTG prediction tool
    f_in = open(LTG_file, "r")
    # is the gene >70% identical to it's blast hit?
    # if so, maybe a contamination?
    gene_to_comment_dict = dict()
    HGT_predicted_gene_set = set([])
    gene_to_HGTspeces_discription_dict = dict()
    gene_to_AI = dict()
    for line in f_in:
        if line.startswith("#"):
            continue
        gene = line.split("\t")[0]
        #call function to format gene name
        gene = split_gff_gene_names(gene)
        comment = line.split("\t")[-1]
        species = line.split("\t")[-3]
        description = line.split("\t")[-2]
        AI = line.split("\t")[7]
        kingdom = line.split("\t")[5]
        gene = gene.replace("ID=", "").split("gene=g")[0]
        gene = gene.rstrip()
        if not ".t" in gene:
            gene = gene.replace("t", ".t")
        HGT_predicted_gene_set.add(gene)
        data_out_formatted2 = "%s\t%s" %(species, description)
        data_out_formatted = "%s\t%s" %(comment.rstrip("\n"),kingdom)
        gene_to_HGTspeces_discription_dict[gene] = data_out_formatted2
        gene_to_comment_dict[gene] = data_out_formatted
        gene_to_AI[gene] = AI
    #print (HGT_predicted_gene_set)
    f_in.close()
    return gene_to_AI

def index_gene_matrix(matrix):
    """function to assign the gene, AI, start stop info
    to a dict.
    The input matrix is: gene AI scaffold start stop
    GROS_g13931	460.5170186	GROS_02504	78	1211
    """
    scaffolf_to_gene_matrix = dict()
    gene_to_scaffold_matrix = dict()
    with open(matrix, "r") as handle:
        #print ("i am here")
        for line in handle:
            if not line.strip():
                continue #if the last line is blank
            if line.startswith("#"):
                continue 
            gene, AI, scaffold, start, stop = line.rstrip("\n").split()
            scaffolf_to_gene_matrix[scaffold]= line
            gene_to_scaffold_matrix[gene]= line.split("\t")
    #return a dictionary. Key[transcript_name], vals are a list containing
    #returns {GROS_02504:GROS_g13931	460.5170186	GROS_02504	78	1211]
    return scaffolf_to_gene_matrix, gene_to_scaffold_matrix
    
def parse_gene_to_scaffold_matrix_vals(vals):
    """function to parse the element of disctionary
    return integer datatype is required
    """
    gene = vals[0]
    AI = vals[1]
    scaffold = vals[2]
    start = int(vals[3])
    stop = int(vals[4])
    return gene, AI, scaffold, start, stop

def parseTransposon_info(transposon_info):
    """function to parse the element of trasposon_file
    . The GFF output from the transposon finding pipeline
    is used here to define the coordinate and annot of the
    transposon. NOTE. The scaffold names may have been changed
    so they fit into the program prior to this. If this brakes.
    check them
    INPUT e.g.
    Mc2825	Repeatmasker-OneCode	transposable_element	26019	26484	2108	-	.	ID=element156452;Name=DNA/hAT-Tip100;Note=?
    Mc2825	Repeatmasker-OneCode	transposable_element	27906	27982	300	-	.	ID=element156453;Name=DNA/RC;Note=?
    Mc2825	Repeatmasker-OneCode	transposable_element	6290	6427	488	-	.	ID=element156450;Name=DNA/hAT;Note=?
    """
    #print transposon_info
    trans_scaffold = transposon_info[0]
    transposon_type = transposon_info[1]
    trans_start = int(transposon_info[3])
    trans_stop = int(transposon_info[4])
    tran_annotation = transposon_info[8]
    tran_annotation = tran_annotation.split(";Note")[0]
    # return the correct datatype
    return trans_scaffold, transposon_type, trans_start, trans_stop, tran_annotation


def disctance_dict(gene_3_primt_dict, gene_5_primt_dict, gene_distances_dict, gene, \
                gene_start, gene_stop, trans_start, trans_stop, AI, \
                   transposon_type, tran_annotation):
    """ function to get three and five prime distances,
    and put info in dictionary for accessing later.
    The function compare the current transpaon distance against values
    already genetared to find the closest transposon."""
    #5 prime distance logic pattern
    if gene_start > trans_stop:
        five_disctance = gene_start - trans_stop
        try:
            # test to see if this gene already has a "closer" 5 prime transposon
            # if new one is closer, set this to the new gene value
            if gene_5_primt_dict[gene] > five_disctance:
                gene_5_primt_dict[gene] = five_disctance
                data_formatted = "five_prime\t%s\t%d\t%s" %(AI,five_disctance,tran_annotation)
                new_name = gene+"_five_prime"
                gene_distances_dict[new_name] = data_formatted
        
        except:
            # we have no data yet for this gene. So add this to the
            # dictionary.
            gene_5_primt_dict[gene] = five_disctance
            data_formatted = "five_prime\t%s\t%d\t%s" %(AI,five_disctance,tran_annotation)
            new_name = gene+"_five_prime"
            gene_distances_dict[new_name] = data_formatted

    #3prime distance logic
    if gene_stop < trans_start:
        three_disctance = trans_start-gene_stop

        try:
            # test to see if this gene already has a "closer" 3 prime transposon
            # if new one is closer, set this to the new gene value

            if gene_3_primt_dict[gene] > three_disctance:
                gene_3_primt_dict[gene] = three_disctance
                data_formatted = "three_prime\t%s\t%d\t%s" %(AI,three_disctance,tran_annotation)
                new_name = gene+"_three_prime"
                gene_distances_dict[new_name] = data_formatted
        except:
            # we have no data yet for this gene. So add this to the
            # dictionary. 
            gene_3_primt_dict[gene] = three_disctance
            data_formatted = "three_prime\t%s\t%d\t%s" %(AI,three_disctance,tran_annotation)
            new_name = gene+"_three_prime"
            gene_distances_dict[new_name] = data_formatted
    # retunr the 3 and 5 prime gene to distances dictionary
    return gene_3_primt_dict, gene_5_primt_dict, gene_distances_dict


def find_nearest_transposon(matrix, LTG_file, gene_gff, transposon_file_in):
    """function to find the nearest transposon"""
    #call the function to index the matrix file
    #note: gene_to_scaffold_matrix is the whole line of the matrix file
    scaffolf_to_gene_matrix, gene_to_scaffold_matrix = index_gene_matrix(matrix)
    #open the transposon GFF
    transposon_file = open(transposon_file_in, "r")
    transposon_read = transposon_file.readlines()
    if matrix:
        matrix_file = open(matrix, "r")
        matrix_file_read = matrix_file.readlines()
    else:
        gene_names, scaffold_to_gene_dict, gene_start_stop_dict =  parse_gff(gene_gff, "r")       

    # empty dictionaries for now
    gene_distances_dict = dict()
    gene_3_primt_dict = dict()
    gene_5_primt_dict = dict()

    for gene_line in matrix_file_read:
        #print gene_line
        if gene_line.startswith("#"):
            continue
        if not gene_line.strip():
            continue #if the last line is blank
        # first element of the files
        gene = gene_line.rstrip().split("\t")[0]
        #note: gene_to_scaffold_matrix is the whole line of the matrix file
        vals = gene_to_scaffold_matrix[gene]
        gene, AI, gene_scaffold, gene_start, gene_stop = parse_gene_to_scaffold_matrix_vals(vals)
        for line in transposon_read:
            if line.startswith("#"):
                continue
            if not line.strip():
                continue #if the last line is blank
            transposon_info = line.rstrip("\n").split("\t")

            #print "gene_line =", gene_line, "transpo = ,", line

            trans_scaffold, transposon_type, trans_start, trans_stop, tran_annotation = parseTransposon_info(transposon_info)
            
            if trans_scaffold == gene_scaffold:
                gene_3_primt_dict, gene_5_primt_dict, gene_distances_dict = disctance_dict(gene_3_primt_dict, \
                                                    gene_5_primt_dict, gene_distances_dict, gene,\
                                                    gene_start, gene_stop, trans_start, trans_stop, AI, \
                                                    transposon_type, tran_annotation)
    return gene_distances_dict
    
                               
################################################################################################################################################

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.3"
    sys.exit(0)


usage = """Use as follows:

$ python transposon_distance.py -t transposon.gff -g gene.gff -o outfile.out

          can use -m instead of -g

script to find the nearest transposon to the genes in both 5 prime and
3prime direction.
This does not take into consideration gene direction
"""

parser = OptionParser(usage=usage)

parser.add_option("-t", "--transposons", dest="transposons_gff", default=None,
                  help="gff file generated for the trapsosons finding",
                  metavar="FILE")
parser.add_option("-g", "--gene", dest="gene_gff", default=False,
                  help="gff file for predicted genes",
                  metavar="FILE")
parser.add_option("-m", "--matrix", dest="matrix", default=False,
                  help="matrix for gene, AI, scaffold, start stop. "
                  " can be used instead of  GFF",
                  metavar="FILE")
parser.add_option("-l", "--ltg", dest="LTG_file", default="LTG_results.out_Alien_index.out",
                  help="LTG_file output from LTG prediction program. "
                  " The file that contains all the AI values "
                  " LTG_results.out_Alien_index.out ",
                  metavar="FILE")
parser.add_option("-o", "--out", dest="outfile", default="transposon_distances_list.out",
                  help="Output filename - default= transposon_distances_list.out",
                  metavar="FILE")
              

(options, args) = parser.parse_args()


#-o
outfile= options.outfile
#-l
LTG_file=options.LTG_file
#-t
transposons_gff= options.transposons_gff
#-g
gene_gff= options.gene_gff
#-m
matrix= options.matrix

##############################################################################################################################
        
# empty string to add the results to
overall_result = ""

# test case
#gene_distances_dict= find_nearest_transposon("./test/test_matrix", "./test/test_trasposon")
# GROS example
gene_distances_dict = find_nearest_transposon("AI_GT_zero_matrix.matrix", "nGr.v1.1_TEs_Thu_Dec__3_09_30_47_2015_Seb.gff")

gene_distances_dict = find_nearest_transposon(matrix, transposons_gff)


for key, vals in gene_distances_dict.items():
    #print vals.split("\t")
    prime, AI,disctance,tran_annotation= vals.split("\t")
    try:
        name = key.split("_three")[0]
    except:
        name = key.split("_five")
    dataformatted = "%s\t%s\t%s\t%s\t%s\n" %(name, prime, AI,disctance,tran_annotation.split(";Name=")[1])
    overall_result = overall_result + dataformatted


overall_list = []

overall_result = overall_result.split("\n")

for i in overall_result:
    overall_list.append(i)


overall_list = sorted(overall_list)

f_out = open("HGT_transposon_distances_list.out", "w")
print >> f_out, "#name\tprime\tAI\tdisctance\ttransposon_annotation"

for i in overall_list:
    #print i
    print >> f_out, i



