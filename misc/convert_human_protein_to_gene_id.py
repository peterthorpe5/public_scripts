""" Title: Convert transcript or gene to gene name
"""

#imports

import os
from sys import stdin,argv
import sys
from collections import defaultdict



def test_line(line):
    """returns true lines. Not comments or blank line"""
    if not line.strip():
        return False  # if the last line is blank
    if line.startswith("#"):
        return False  # comment line
    return line


def parse_names():
    """we need o convert the trans or gene name to the func name
    >ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|
    OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|

    reruns to default dict e.g.
    T  = transcript
    G = gene
    [ENST] = DDX11L1
    [ENSG] = DDX11L1

    """
    trans_to_gene = defaultdict(str)
    gene_to_trans = defaultdict(list)
    tran_to_func = defaultdict(str)
    gene_to_func = defaultdict(str)
    #
    trans_to_gene_id = defaultdict(str)
    gene_to_gene_id = defaultdict(str)
    gene_id_to_func = defaultdict(str)
    f_in = open("headers", "r")
    for line in f_in:
        if test_line(line):
            data = line.split("|")
            transcript = data[0]
            transcript = transcript.replace(">", "")
            gene = data[1]
            gene_id = data[5]
            func = data[7]
            # populate the dictionaries
            trans_to_gene[transcript] = gene
            gene_to_trans[gene].append(transcript)
            tran_to_func[transcript] = func
            gene_to_func[gene] = func
            #            
            trans_to_gene_id[transcript] = gene_id
            gene_to_gene_id[gene] = gene_id
            gene_id_to_func[gene_id] = func
    return trans_to_gene, gene_to_trans, tran_to_func, \
           gene_to_func, trans_to_gene_id, \
           gene_to_gene_id, gene_id_to_func
            

def parse_file(infile, outfile, trans_to_gene,
               gene_to_trans, tran_to_func,
               gene_to_func, trans_to_gene_id, 
               gene_to_gene_id, gene_id_to_func):
    f_in = open(infile, "r")
    f_out = open(outfile, "w")

    for line in f_in:
        if test_line(line):
            if line.startswith("sampleA"):
                f_out.write("\t" + "gene_id" + "\t" + line)
                continue
            data = line.split("\t")
            subject = data[0]
            if "transcript" in infile:
                gene_id = trans_to_gene_id[subject]
                func = tran_to_func[subject]
                out_data = "%s\t%s %s " % (subject, gene_id, func)
                line = line.replace(subject, out_data)
                f_out.write(line)
                
            if "gene" in infile:
                gene_id = gene_to_gene_id[subject]
                func = gene_to_func[subject]
                out_data = "%s %s %s " % (subject, gene_id, func)
                line = line.replace(subject, out_data)
                f_out.write(line)
    f_in.close()
    f_out.close()
    

if __name__ == '__main__':
    trans_to_gene, gene_to_trans, tran_to_func, \
           gene_to_func, trans_to_gene_id, \
           gene_to_gene_id, gene_id_to_func = parse_names()
    # call the function to get a list of results wanted
    directory = "."
    # get all folder names os.walk(directory)
    for dirpath, subdirs, files in os.walk(directory):
        for x in files:
    
            if x.endswith(".pdf"): continue
            if x.endswith("subset"):
                wanted_file = (os.path.join(dirpath, x))
                outfile = wanted_file + "_RENAMED"
                #print(wanted_file)
                parse_file(wanted_file, outfile, trans_to_gene,
                           gene_to_trans, tran_to_func,
                           gene_to_func, trans_to_gene_id, 
                           gene_to_gene_id, gene_id_to_func)
