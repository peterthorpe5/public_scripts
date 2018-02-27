#!/usr/bin/env python

########################################################################################################
#####: TITLE: Get the top blast hit from tab file. Get Kingdom hits. Taxonmy filter  ###############


"""
script to return  the top blast hit from tab file (via two methods: One the assuming order in the blast output
and 2) explicitly looking for the hit with the greatest bit score). The distribution of the kingdom
for the top hits is also returned in a file. 

#to do taxonomy filter. The code is there, just need to implement this as an option in the main function 
"""
#author: Peter Thorpe September 2015. The James Hutton Insitute, Dundee, UK.

#biopython imports

import time
#os imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser


##################################################################################################
def wanted_genes(blast_file):
    "function to retunr a list of wanted genes from file"
    wanted = open(blast_file, "r")
    names = wanted.readlines()
    blast_data = [line.rstrip() for line in names
              if line.strip() != "" if not line.startswith("#")]
    wanted.close()
    #print "wanted_data :", blast_data
    return blast_data

def is_number(s):
    "check cloumn if a float number"
    try:
        float(s)
        return s-1
    except ValueError:
        print "coloumn may not be the correct bit score coloumn. Please check"
        return False
    
###########################################################################################################################################################################################

  
def get_top_blast_hit(wanted_data, bit_score_column, outfile):
    "parse through file and get top hit. Prints to a file reduced info."
    blast_data = wanted_genes(in_file)
    bit_score_column = int(bit_score_column)-1

    got_list = set([])
    outfile_name = outfile+"_based_on_order_tax_king.tab"
    f = open(outfile_name, "w")
    for line in blast_data:
        name = line.split("\t")[0]
        #print name
        description = line.split("\t")[bit_score_column+1]
        line = line.rstrip("\n")
        if not name in got_list:
            wanted_info = "%s\n" %(line)
            f.write(wanted_info)
            got_list.add(name)
    f.close()


###########################################################################################################################################################################################


def parse_NCBI_nodes_tab_file(folder):
    """this is a function to open nodes.dmp from the NCBI taxonomy
database and find the parent child relationship....returns a
disctionary for later use"""

    #open file - read.
    #nodes.dmp - this file is separated by \t|\t
    #empty dictionary to add to parent and child (keys,vals) to
    tax_dictionary = {}

    #nodes.dmp files goes: child, parent, etc
    #merged.dmp file goes: old, new
    #In both cases, can take key as column 0 and value as column 1
    for filename in ["nodes.dmp", "merged.dmp"]:
        with open(os.path.join(folder, filename)) as handle:
            for line in handle:
                tax_info = line.replace("\n", "\t").split("\t|\t")
                #first element
                parent = tax_info[1]
                #second element
                child = tax_info[0]
                #add these to the dictionary {parent:child}
                tax_dictionary[child]= parent
    #print tax_dictionary    
    return tax_dictionary


###########################################################################################################################################################################################

def test_if_id_is_nematoda(tax_dictionary, tax_id_of_interst,final_tx_id_to_identify_up_to, tax_to_filter_out):
    """function to get a list of tax id of interest from the tax_dictionary
    which is produced in the parse_function (parse_NCBI_nodes_tab_file)
    nodes.dmp file. . The tax id
    are added to a list for later use

    """
    tax_id_of_interst= tax_id_of_interst.strip()
    if tax_id_of_interst == "0":
        raise ValueError("0 is an unknown ID, going to assing 32644 to it instead")
        tax_id_of_interst ="32644" #assign an unknown tax_id
    if tax_id_of_interst == "N/A":
        raise ValueError("N/A as taxonomy ID")
    #get the "master" parent id
    parent = tax_dictionary[tax_id_of_interst]
    #print parent

    while True:
        #print "parent = ", parent, "\n"
        parent = tax_dictionary[parent]
        if tax_id_of_interst == "N/A":
            raise ValueError("N/A as taxonomy ID")
        #32630 is a synthetic organism
        if parent == "32630":#32630
            return "In_filter_out_tax_id"
            break            
        if parent == tax_to_filter_out:
            return "In_filter_out_tax_id"
            break
        if parent == final_tx_id_to_identify_up_to:
            #print "......................... im here"
            return True
        elif parent == "1":
            # Reached the root of the tree
            return False



def get_genus_count(genus_dict, blast_line, sci_name_column="15"):
    """this function count the distribution of the genus for the top hit"""
    sci_name_column = int(sci_name_column)-1
    scinetific_name = blast_line[sci_name_column]
    try:
        genus = scinetific_name.split()[0]
    except:
        return genus_dict
    try:
        genus_dict[genus]+=1
    except:
        genus_dict[genus]=1
    return genus_dict

##############################################################################################
def parse_blast_tab_file(in_file, bit_score_column, outfile):
    """this is a function to open up a tab file blast results, and
    produce the percentage of kingdom blast hit based on the top
    blast match"""

    blast_data = wanted_genes(in_file)
    get_top_blast_hit(blast_data, bit_score_column, outfile)
    #open files, read and write.
    blast_file = open (in_file, "r")
    out_file = open(outfile,"w")
    bit_score_column = int(bit_score_column)-1

    
    #set of blast_file_entry gene names
    blast_file_entry_Genes_names = set([])
    kingdoms = set("")

    # dictionary of all the kingdoms in our blast file
    kingdoms_handles_counts = {'Eukaryota':0, 'N/A':0, 'Bacteria;Eukaryota':0, \
                               'Archaea;Eukaryota':0, 'Virus':0, 'Bacteria;Viruses':0,\
                               'Eukaryota;Viruses':0, 'Archaea':0, 'Bacteria':0,\
                               'Unclassified':0}
    
    #this is out list of so called top matches which we will append and remove as applicable
    top_hits = []
    #current bit score value "to beat"
    current_bit_score = float(0.0)
    last_gene_name = ""
    last_blast_line = ""
    
    for line in blast_file:
        if line.startswith("#"):
            continue
        if not line.strip():
            continue
        #print line
        blast_line = line.rstrip("\n").split("\t")
        #names of the query seq
        blast_file_entry_Genes = blast_line[0]
        #print blast_file_entry_Genes
        bit_score = float(blast_line[bit_score_column])
        kings_names = blast_line[-2]
        
        ##############################################################################
        #first block: if the names are the same, is the new bit score more?
        if blast_file_entry_Genes == last_gene_name:
            #print "im here"
            if bit_score > current_bit_score:
                #print "current_bit_score", current_bit_score
                current_bit_score = bit_score
                #print "current_bit_score", current_bit_score
                #remove the last entry if so and put the new one in
                del top_hits[-1]
                top_hits.append(blast_line)

                
        #############################################################################
        # second block: if the name is new, put it in the name set.
        # use this bit score as the new one to "beat"

        #print current_bit_score
        if not blast_file_entry_Genes in blast_file_entry_Genes_names:
            #print ".......should be first line"
            blast_file_entry_Genes_names.add(blast_file_entry_Genes)
            current_bit_score = bit_score
            top_hits.append(blast_line)

        ############################################################################
        # assign value to the variables for testing in the new batch of for loops
        
        last_gene_name = blast_file_entry_Genes
        last_blast_line = line
        
    genus_dict = dict()

    total_blast_hit_count = 0 
    for i in top_hits:
        genus_dict = get_genus_count(genus_dict, i)
        total_blast_hit_count = total_blast_hit_count+1
        king_name = i[-2]
        print "KING", king_name
        kingdoms_handles_counts[king_name]+=1
        new_line = ""
        for element in i:
            new_line = new_line+element+"\t"
            
        data_formatted = new_line.rstrip("\t")+"\n"
        out_file.write(data_formatted)



    #for blast_file_entry_Genes, bit_score, kings_names in top_hits:
        #print >> out_file, "%s\t%s\t%s" % (blast_file_entry_Genes, bit_score, kings_names)
        #kingdoms_handles_counts[kings_names]+=1

    print "Kingdom hit distribution of top hits = ", kingdoms_handles_counts
    print "number with blast hits =", total_blast_hit_count
    #print "genus distirbution =", genus_dict

    top_hits_out_king = open("kingdom_top_hits.out", "w")
    file_tile = "#top kingdom hit for %s\n" %(in_file)
    top_hits_out_king.write(file_tile)

    top_hits_out_genus = open("Genus_distribution_top_hits.out", "w")
    file_tile = "#Genus_of_top hits for %s\n" %(in_file)
    top_hits_out_genus.write(file_tile)

    for kingdom, count in kingdoms_handles_counts.items():
        data_form = "%s:\t%i\n" %(kingdom, count)
        top_hits_out_king.write(data_form)

    for genus, count in genus_dict.items():
        data_form = "%s:\t%i\n" %(genus, count)
        top_hits_out_genus.write(data_form)

    top_hits_out_king.close()
    out_file.close()
    top_hits_out_genus.close()
    
    return kingdoms_handles_counts




#####################################################################################################
start_time=time.time()
###################################################################################################


if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.2"
    sys.exit(0)


usage = """Use as follows:

$ python top_BLAST_hit.py -i in.tab -b bit_score coloumn -o out_file 

This iterates through a tabular blast output
and returns the top blast hit based on the greatest bitscore (-b colomn in file - default is 12 12-1 computer counting)

Usually the order of blast hit is preservedm, this script explicitly looks for the best, checking the order has not been altered. 




script to return  the top blast hit from tab file (via two methods: One the assuming order in the blast output
and 2) explicitly looking for the hit with the greatest bit score). The distribution of the kingdom
for the top hits is also returned in a file. 

#to do taxonomy filter. The code is there, just need to implement this as an option in the main function 

"""

parser = OptionParser(usage=usage)


parser.add_option("-i", "", dest="in_file", default=None,
                  help="in file")
parser.add_option("-b", "--bit_score", dest="bit_score_column", default="12",
                  help="bit_score_column")

parser.add_option("-o", "--output", dest="out_file", default=None,
                  help="Output filename",
                  metavar="FILE")


(options, args) = parser.parse_args()

in_file = options.in_file
bit_score_column = options.bit_score_column
outfile = options.out_file


#if len(args) < 1:
    #stop_err("Expects no argument, one input filename")



#get_top_blast_hit(in_file, bit_score_column, outfile)

parse_blast_tab_file(in_file, bit_score_column, outfile)


end_time=time.time()
#print 'that took, %.3f' %(end_time - start_time)

#print 'done'
