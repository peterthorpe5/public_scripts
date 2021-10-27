# script to return gene counts for GOI
# author Pete Thorpe 2021 Oct
# imports

import os
import sys
from optparse import OptionParser
import datetime
import logging
import logging.handlers
import time
from collections import defaultdict


def test_line(line):
    """returns true lines. Not comments or blank line"""
    if not line.strip():
        return False  # if the last line is blank
    if line.startswith("#"):
        return False  # comment line
    return line


def get_trans_to_gene(gene_to_trans, logger):
    """func take in a file with:.
    # GENE  gene_symbol trnacripts
    ENSMUSG00000092210.1	A930009A15Rik	ENSMUST00000173620.1
    returns a dictionaries"""
    f_in = open(gene_to_trans, "r")
    gene_to_transcript = defaultdict(list)
    transcript_to_gene = defaultdict(str)
    gene_to_symbol = defaultdict(str)
    for line in f_in:
        if test_line(line):
            data = line.split()
            gene = data[0]
            symbol = data[1]
            gene_to_symbol[gene] = symbol
            short_gene_name = gene.split(".")[0]
            # the result file look like this ENSMUSG00000029273
            # rather than ENSMUSG00000029273.13
            gene_to_symbol[short_gene_name] = symbol
            for entry in data:
                if "ENSMUST" in entry:
                    transcript = entry.rstrip()
                    gene_to_transcript[gene].append(transcript)
                    short_gene_name = gene.split(".")[0]
                    gene_to_transcript[short_gene_name].append(transcript)
                    
                    transcript_to_gene[transcript] = gene
                    transcript_to_gene[transcript] = short_gene_name
    return gene_to_transcript, transcript_to_gene, gene_to_symbol
                    
                    
            
def parse_genes(genes, logger):
    """func take in a file with a list og gene names. 
    returns a set"""
    f_in = open(genes, "r")
    input_count = 0
    goi_set = set([])
    for line in f_in:
        input_count = input_count + 1
        if test_line(line):
            line = line.strip()
            line = line.upper()
            goi_set.add(line)
    info = "input GOI list was %d genes" % input_count
    logger.info(info)
    info = "non-redundant GOI list was %d genes" % len(goi_set)
    logger.info(info)   
    return(goi_set)



def parse_DE(gene_to_transcript, transcript_to_gene, 
             gene_to_symbol, transcript, gene_set,
             outfile, logger):
    """fucn take in the DE results file and compare to
    three dictionaies already create. """
    outfmt = "looking at %s" % transcript
    logger.info(outfmt)
    trans_set = set([])

    f_in = open(transcript, "r")
    f_out = open(outfile, "w")
    DE_transcript_count = 0
    logger.info("DE genes")
    logger.info(gene_set)
    for line in f_in:
        if test_line(line):
            if line.startswith("Row.names"):
                continue
            DE_transcript_count = DE_transcript_count + 1
            data = line.split()
            name = data[0]
            trans_set.add(name)
            gene = transcript_to_gene[name]
            symbol = gene_to_symbol[gene]
            # get the corresponding gene name
            if gene in gene_set:
                outfmt = ("DE transcript \t%s\thas DE gene \t%s\tSYMBOL=\t%s\tSUCCESS" % (name, gene, symbol))
                logger.info(outfmt)
            if gene not in gene_set:
                outfmt = ("DE transcript\t%s\tdoes NOT have DE gene \t%s\tSYMBOL=\t%s" % (name, gene, symbol))
                logger.info(outfmt)
    outfmt = "DE transcript found =  %d" % DE_transcript_count
    logger.info("DE transcript")
    logger.info(trans_set)
    logger.info(outfmt)
                
                

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

$ python filter....py -h

"""

parser = OptionParser(usage=usage)

parser.add_option("-g", dest="gene",
                  default="DE_gene_veh_vs_ALDO",
                  help="list of DE genes")

parser.add_option("-t", dest="transcript",
                  default="aldo_vs_veh.GLM.passage_as_batch_edgeR.DE_LOGFC_1.0_FDR_0.05",
                  help="De transcript file",
                  metavar="FILE")

parser.add_option("-m", dest="map",
                  default="mus_genes_to_transcripts.txt",
                  help="gene to trascript map",
                  metavar="FILE")
                  
parser.add_option("-o", dest="outfile",
                  default="test.out",
                  help="out file for results")



(options, args) = parser.parse_args()

gene_to_trans = options.map
gene = options.gene
transcript = options.transcript

outfile = options.outfile

logfile = "%s_COMPARED_TO_%s_RESULTS.log" % (gene, transcript)
# Run as script
if __name__ == '__main__':
    start_time = time.time()
    # Set up logging
    logger = logging.getLogger('compare...py: %s'
                               % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)
    try:
        logstream = open(logfile, 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        logger.error("Could not open %s for logging" %
                     logfile)
        sys.exit(1)
    # get gene of interests
    gene_to_transcript, transcript_to_gene, \
                gene_to_symbol = get_trans_to_gene(gene_to_trans, logger)
    # get the DE genes
    gene_set = parse_genes(gene, logger)
    # parse the DE results 
    parse_DE(gene_to_transcript, transcript_to_gene, 
             gene_to_symbol, transcript, gene_set,
             outfile, logger)
            
