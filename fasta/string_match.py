#!/usr/bin/env python
# coding: utf-8
# script to search for perefect peptide matches in a fasta db
# author Pete Thorpe 2022 April
# imports


#from fuzzysearch import find_near_matches
from collections import defaultdict
import collections
#import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from optparse import OptionParser
import datetime
import logging
import logging.handlers
import time
import sys
import re



def get_peps(infile, logger):
   """func to get peps: returns a set
   infile is csv. Peptide coloumn is called Peptide"""
   data = pd.read_csv(infile) #, skiprows = 1)
   peptides = list(data["Peptide"])
   d = ("we have an original list of %d peptides" % (len(peptides)))
   logger.info(d)

   pep_set = set([])
   for pep in peptides:
      pep = pep.split("(")[0] # there is extra weirdness on some.
      pep_set.add(pep)
   d =("we now have a non-reducndant list of %d peptides" % (len(pep_set)))
   logger.info(d)
   return pep_set


def seq_getter(fasta, pep_set, logger):
    """this is a function to open up a fasta file and
    returns those which are
    contained within the set.
    returns a results dict: [pep] = list of hist
    list is coordinate. seq.id"""
    reults_dict = defaultdict(list)

    for seq_record in SeqIO.parse(fasta, "fasta"):
         for pep in pep_set:
            if pep in str(seq_record.seq):
               #print("yes", seq_record.id)
               # print(pep, str(seq_record.seq)[:10])
               result = str(seq_record.seq).find(pep)
               out_result = "\t".join([str(result), seq_record.id])
               reults_dict[pep].append(out_result)
    return reults_dict



if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

$ python ....py --fasta db.fasta --csv peptide.csv

"""

parser = OptionParser(usage=usage)

parser.add_option("--fasta", dest="fasta",
                  default="DB_SARS2_6frame_trans_Novel_ORF_SARSAA_Uniprot_WITH_ORF_non_redundant.fasta",
                  help="db fasta file")


parser.add_option("--csv", dest="csv",
                  default=None,
                  help="csv file",
                  metavar="FILE")

parser.add_option("--out", dest="out",
                  default="results.txt",
                  help="outfile for the results",
                  metavar="FILE")



(options, args) = parser.parse_args()

fasta = options.fasta
csv = options.csv
outfile = options.out



logfile = "peptide_searching.log"
# Run as script
if __name__ == '__main__':
    start_time = time.time()
    # Set up logging
    logger = logging.getLogger('string_match.py: %s'
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
    pep_set = get_peps(csv, logger)
    results_dict = seq_getter(fasta, pep_set, logger)
    
    f = open(outfile, 'w')
    title = "#Peptide\tcoordinate_start\tSeq_id\tcoordinate_start2\tSeq_id\n"
    f.write(title)
    for pep, hits in results_dict.items():
       out = ""
       
       for i in hits:
          out = out + i + "\t"
       out = pep + "\t" + out + "\n"
       f.write(out)
    f.close()


     
          


         

