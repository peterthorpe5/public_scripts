"Ttitle: script to search for self generated list of gene, to their nearest transposon."
#why: The question is: Are HGT genes closeer to transposons?
# need to test the real value against a random set to see if statistically different.
###############################################


#author: Peter Thorpe September 2015. The James Hutton Insitute, Dundee, UK.

#imports
import random
import os
import sys
from optparse import OptionParser
import datetime


##############################################
#dictionary for trans assignemnt:

transposon_to_class_dict = {'LTR/ERVK':'LTR/ERVK','LTR/Roo':'LTR/Roo', 'SINE/ID':'SINE/ID','DNA/TcMar-Tc1': '2', 'DNA/TcMar-Tc2': '18',\
                            'LINE/CR1': '33', 'DNA/TcMar-Tc4': '17', 'LTR/Pao': '22', 'DNA/Chapaev': '26', 'LINE/L2': '31', 'LINE/L1': '16',\
                            'LINE/Unclassified': '14', 'LINE/R1': '35', 'DNA/hAT-Tag1': '25', 'LINE/R2-NeSL': '38', 'DNA/hAT-hATx': '23',\
                            'SINE?': '7', 'DNA/hAT': '10', 'DNA/TcMar-Mariner': '20', 'DNA/Transib': '8',\
                            'LINE/RTE-RTE': '34', 'DNA/Maverick': '21', 'DNA/Helitron': '19', 'DNA/Ginger': '5',\
                            'DNA/Sola': '37', 'SINE/Alu': '29', 'LINE/Jockey': '30', 'DNA/Merlin': '24', 'LTR/Gypsy': '3',\
                            'DNA/hAT-Charlie': '15', 'DNA/hAT-Ac': '1', 'DNA/ISC1316': '36', 'DNA/EnSpm': '4',\
                            'DNA/TcMar-Tigger': '27', 'DNA/Unclassified': '9', 'LTR/Unclassified': '11', 'DNA/hAT-Blackjack': '28',\
                            'DNA/TcMar-Pogo': '12', 'LTR/Copia': '6', 'SINE/5S': '32', 'DNA/MuLE-MuDR': '13'}

transposon_count_dict = {'LTR/ERVK':0, 'LTR/Roo': 0,'SINE/ID': 0, 'DNA/TcMar-Tc1': 0, 'DNA/TcMar-Tc2': 0, 'LINE/CR1': 0,\
                         'DNA/TcMar-Tc4': 0, 'LTR/Pao': 0, 'DNA/Chapaev': 0, 'LINE/L2': 0, 'LINE/L1': 0,\
                         'LINE/Unclassified': 0, 'LINE/R1': 0, 'DNA/hAT-Tag1': 0, 'LINE/R2-NeSL': 0,\
                         'DNA/hAT-hATx': 0, 'SINE?': 0, 'DNA/hAT': 0, 'DNA/TcMar-Mariner': 0, 'DNA/Transib': 0,\
                         'LINE/RTE-RTE': 0, 'DNA/Maverick': 0, 'DNA/Helitron': 0, 'DNA/Ginger': 0, 'DNA/Sola': 0,\
                         'SINE/Alu': 0, 'LINE/Jockey': 0, 'DNA/Merlin': 0, 'LTR/Gypsy': 0, 'DNA/hAT-Charlie': 0,\
                         'DNA/hAT-Ac': 0, 'DNA/ISC1316': 0, 'DNA/EnSpm': 0, 'DNA/TcMar-Tigger': 0, 'DNA/Unclassified': 0,\
                         'LTR/Unclassified': 0, 'DNA/hAT-Blackjack': 0, 'DNA/TcMar-Pogo': 0, 'LTR/Copia': 0, 'SINE/5S': 0,\
                         'DNA/MuLE-MuDR': 0}

#####################################################################################################################################################################

def random_name_generator(name_file):
    "function to get a random list of names"
    "name_file: just a list of all gene names. cut this from the GFF"
    f_in = open(name_file, "r")
    wanted_data = [line.replace("\t", "").rstrip("\n") for line in f_in
              if line.strip() != ""]
    name_set = set([])
    for i in range(0,519):
        gene = random.choice(wanted_data)
        if not gene.startswith("#"):
            name_set.add(gene)
    #print wanted_data
    #returns a set of genes
    return name_set

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


def find_nearest_transposon(matrix, transposon_file_in):
    """function to find the nearest transposon"""
    #call the function to index the matrix file
    #note: gene_to_scaffold_matrix is the whole line of the matrix file
    scaffolf_to_gene_matrix, gene_to_scaffold_matrix = index_gene_matrix(matrix)
    #open the transposon GFF
    transposon_file = open(transposon_file_in, "r")
    transposon_read = transposon_file.readlines()
    matrix_file = open(matrix, "r")
    matrix_file_read = matrix_file.readlines()

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
##############################################################################################################################################################################    
                               

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.3"
    sys.exit(0)


usage = """Use as follows:

$ python transposon_distance.py -t transposon.gff -g gene.gff -i 100 -o outfile.out

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
parser.add_option("--gene_list", dest="gene_list", default=False,
                  help="list of gene names to randomly select from",
                  metavar="FILE")
parser.add_option("-m", "--matrix", dest="matrix", default=False,
                  help="matrix for gene, AI, scaffold, start stop. "
                  " can be used instead of  GFF",
                  metavar="FILE")
parser.add_option("-i", "--iterations", dest="iterations", default=False,
                  help="number of random iterations to perform "
                  " when generateing random lists",
                  metavar="FILE")
parser.add_option("-o", "--out", dest="outfile", default="transposon_distances_list.out",
                  help="Output filename - default= transposon_distances_list.out",
                  metavar="FILE")
              

(options, args) = parser.parse_args()


#-o
outfile= options.outfile
#-t
transposons_gff= options.transposons_gff
#-g
gene_gff= options.gene_gff
#-m
matrix= options.matrix
#-i
iterations=options.iterations
# --gene_list
gene_list = options.gene_list

##############################################################################################################################

                
        
for repeat in range(0,iterations):
    overall_result = ""
    # set up a dict for counting. 
    transposon_count_dict = {'DNA/TcMar-Tc1': 0, 'DNA/TcMar-Tc2': 0, 'LINE/CR1': 0, 'DNA/TcMar-Tc4': 0, \
                             'DNA/TcMar-Ant1': 0, 'LTR/Pao': 0, 'LTR/ERVK': 0, 'DNA/Chapaev': 0, 'LINE/L2': 0, \
                             'LINE/L1': 0, 'LINE/R1': 0, 'DNA/hAT-Tag1': 0, 'LTR/Roo': 0, 'LINE/R2-NeSL': 0,\
                             'SINE/ID': 0, 'DNA/hAT-hATx': 0, 'SINE?': 0, 'DNA/hAT': 0, 'DNA/TcMar-Mariner': 0,\
                             'DNA/Transib': 0, 'LINE/RTE-RTE': 0, 'DNA/P_element': 0, 'DNA/Maverick': 0,\
                             'DNA/Helitron': 0, 'LINE/Unclassified': 0, 'DNA/Sola': 0, 'SINE/Alu': 0,\
                             'LINE/Jockey': 0, 'DNA/Merlin': 0, 'LTR/Gypsy': 0, 'DNA/hAT-Charlie': 0,\
                             'DNA/hAT-Ac': 0, 'DNA/Ginger': 0, 'DNA/ISC1316': 0, 'DNA/EnSpm': 0, 'YR/Ngaro': 0,\
                             'LTR/Copia': 0, 'DNA/Unclassified': 0, 'LTR/Unclassified': 0, 'DNA/hAT-Blackjack': 0, \
                             'DNA/TcMar-Pogo': 0, 'LTR/ERV1': 0, 'DNA/TcMar-Tigger': 0, 'SINE/5S': 0, 'DNA/MuLE-MuDR': 0}
    #gene_distances_dict= find_nearest_transposon("./test/test_matrix", "./test/test_trasposon")
    #gene_distances_dict = find_nearest_transposon("format_for_py_script.out_with_fake_AI.txt", "not_alien.names", \
                                            #"nGr.v1.1_TEs_Thu_Dec__3_09_30_47_2015_Seb.gff")
    gene_distances_dict = find_nearest_transposon(matrix, gene_list, transposons_gff)
    #loop over the dictionary to get the values out. 
    for key, vals in gene_distances_dict.items():
        #print vals.split("\t")
        primt, AI,disctance,tran_annotation= vals.split("\t")
        try:
            name = key.split("_three")[0]
        except:
            name = key.split("_five")
        transposon_type = tran_annotation.split(";Name=")[1]
        #count the number of times we see this transposon
        transposon_count_dict[transposon_type] +=1 
        #format to write out
        dataformatted = "%s\t%s\t%s\t%s\t%s\n" %(name.split("_five_prime")[0], primt, AI,disctance,tran_annotation.split(";Name=")[1])
        overall_result = overall_result + dataformatted

    
    overall_list = []

    overall_result = overall_result.split("\n")

    for i in overall_result:
        overall_list.append(i)


    overall_list = sorted(overall_list)
    out_name = "Random_transposon_distances_list_%d.out" % (repeat)
    
    transposon_file_name = "Random_transposon_classes_%d.out" % (repeat)
    transposon_dictOut = open(transposon_file_name, "w")
    
    print >> transposon_dictOut,"#trans_type\tcount\tclassed"
    for key, vals in transposon_count_dict.items():
        trans_type = key
        count = vals
        try:
            classed = transposon_to_class_dict[trans_type]
        except:
            classed = trans_type
        #print "%s\t%d\t%s" %(trans_type, count, classed)
        print >> transposon_dictOut,"%s\t%d\t%s" %(trans_type, count, classed)

    f_out = open(out_name, "w")
    print >> f_out, "#name\tprime\tAI\tdisctance\ttransposon_annotation"

    for i in overall_list:
        #print i
        print >> f_out, i
    print "done file %d" % (repeat)
    f_out.close()
    transposon_dictOut.close()

print "finished"



