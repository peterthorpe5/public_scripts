""" Title: Convert transcript or gene to gene name
"""

#imports

import os
from sys import stdin,argv
import sys
from collections import defaultdict
from optparse import OptionParser
import gzip
import logging
import logging.handlers
import time


def test_line(line):
    """returns true lines. Not comments or blank line"""
    if not line.strip():
        return False  # if the last line is blank
    #if line.startswith("#"):
        #return False  # comment line
    return line


def parse_names():
    """we need o convert the trans or gene name to the func name
    >lcl|NC_061794.1_cds_XP_047226078.1_1 [gene=LOC124871109]
    [db_xref=GeneID:124871109]
    [protein=uncharacterized protein LOC124871109]
    [protein_id=XP_047226078.1]
    [location=complement(join(4246..4643,4752..
    4970,5054..5283,5499..5696,5774..5910,6005..6121,623
    1..6438,6567..6634))] [gbkey=CDS]


    reruns to default dict e.g.
    T  = transcript
    G = gene
    [XP_047213791.1_28] = phc2a


    """
    trans_to_gene = defaultdict(str)
    gene_to_trans = defaultdict(list)
    tran_to_func = defaultdict(str)
    gene_to_func = defaultdict(str)
    protein_id_to_gene = defaultdict(str)
    #
    trans_to_gene_id = defaultdict(str)
    gene_to_gene_id = defaultdict(str)
    gene_id_to_func = defaultdict(str)
    with gzip.open("headers.gz", mode="rt") as f_in:
        for line in f_in:
            if test_line(line):
                data = line.split("|")
                transcript = data[1]
                transcript_NC = transcript.split("_cds_")[0]
                protein_id = line.split("_cds_")[1]
                protein_id = protein_id.split(" [gene=")[0]
                protein_id = "_".join(protein_id.split("_")[:-1])

                gene = line.split("[gene=")[1]
                gene = gene.split("]")[0]
                gene_id = gene
                func = line.split("[protein=")[1]
                func = func.split("]")[0]
                # populate the dictionaries
                trans_to_gene[transcript] = gene
                gene_to_trans[gene].append(transcript)
                tran_to_func[transcript] = func
                gene_to_func[gene] = func
                protein_id_to_gene[protein_id] = gene
                #            
                trans_to_gene_id[transcript] = gene_id
                gene_to_gene_id[gene] = gene_id
                gene_id_to_func[gene_id] = func
    return trans_to_gene, gene_to_trans, tran_to_func, \
           gene_to_func, trans_to_gene_id, \
           gene_to_gene_id, gene_id_to_func, protein_id_to_gene
            

def parse_file(infile, outfile, trans_to_gene,
               gene_to_trans, tran_to_func,
               gene_to_func, trans_to_gene_id, 
               gene_to_gene_id, gene_id_to_func,
               protein_id_to_gene):
    f_in = open(infile, "r")
    f_out = open(outfile, "w")

    for line in f_in:
        if test_line(line):
            if line.startswith("#query"):
                new_line = line.replace("query", "query\tgene_id")
                f_out.write(new_line)
                continue
            if line.startswith("#"):
                f_out.write(line)
                continue

            data = line.split("\t")
            subject = data[0]
            gene_id = protein_id_to_gene[subject]
            #print(subject, gene_id)
            # func = gene_to_func[gene_id]
            out_data = "%s\t%s" % (subject, gene_id)
            line = line.replace(subject, out_data)
            f_out.write(line)
                
    f_in.close()
    f_out.close()
    

usage = """Use as follows:

$ python convert....py -h 

"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="infile",
                  default=None,
                  help="can fill this in later for sepcfic files")


(options, args) = parser.parse_args()

infile = options.infile



if __name__ == '__main__':
    trans_to_gene, gene_to_trans, tran_to_func, \
           gene_to_func, trans_to_gene_id, \
           gene_to_gene_id, gene_id_to_func, \
           protein_id_to_gene = parse_names()
    # first_three_items = list(protein_id_to_gene.items())[:3]
    # print(first_three_items)

    # call the function to get a list of results wanted
    directory = "."
    # get all folder names os.walk(directory)
    for dirpath, subdirs, files in os.walk(directory):
        for query_file in files:
    
            if query_file.endswith(".pdf"): continue
            if query_file.endswith("annotations"):
                wanted_file = (os.path.join(dirpath, query_file))
                outfile = wanted_file + "_RENAMED"
                #print(wanted_file)
                parse_file(wanted_file, outfile, trans_to_gene,
                           gene_to_trans, tran_to_func,
                           gene_to_func, trans_to_gene_id, 
                           gene_to_gene_id, gene_id_to_func,
                           protein_id_to_gene)
