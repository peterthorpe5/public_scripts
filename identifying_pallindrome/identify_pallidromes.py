#!/usr/bin/env python
#Title: This is a "master" script that import a load of function and
#runs the overall pallidrome identification

#imports
import time
import sys
import os
from optparse import OptionParser
#import functions
from palindrome_library import load_genome_sequence, \
                               iteration_length, \
                               write_palindromes_to_file, \
                               write_shuffled_genome

#######################################################################

usage = """

Usage:
./identify_pallindromes .py --fasta /path/file.fasta --folder (where_ .fasta is, defulat, current working dir) --base_name (species)
    --min_pallindrome_len (default 7) --max_pallindrome_len (default 1000)

This is a master script to call the palindrome library. This will identify
pallindromes the come under a certain probability of occuring at random and
allows missmtaches...

The search will look for pallindromes from min_length to max_length
"""
parser = OptionParser(usage=usage)

parser.add_option("--fasta", dest="fasta", default=None,
                  help="fasta for the pallindrom analysis"
                  " this has to be one contig per file - sorry",
                  metavar="file")
parser.add_option("--folder", dest="folder", default=os.getcwd(),
                  help="folder name where the fasta file is"
                  " results will be put into here",
                  metavar="folder")
parser.add_option("--base_name", dest="base_name", default=None,
                  help="base_name for the species")
parser.add_option("--mask_NNN", dest="mask_NNN", default="False",
                  help="do not return any N containing sequences")
parser.add_option("--repeats", dest="repeats", default=500,
                  help="number of repeats for random interations")
parser.add_option("--min_pallindrome_len", dest="min_pallindrome_len", default=7,
                  help="min length to look for pallindrome. Default 7")
parser.add_option("--max_pallindrome_len", dest="max_pallindrome_len", default=1000,
                  help="max len to look for pallindromes for. Default 1000")


(options, args) = parser.parse_args()

#--fasta
fasta = options.fasta
#--folder
folder = options.folder
#--base_name
base_name = options.base_name
#--max_pallindrome_len
max_pallindrome_len = options.max_pallindrome_len
#--min_pallindrome_len
min_pallindrome_len = options.min_pallindrome_len
#--mask_NNN
mask_NNN = options.mask_NNN
#--repeats
repeats = options.repeats

##########################################################################


dest_dir = os.path.join(folder, 'random')
try:
    os.makedirs(dest_dir)
except OSError:
    print "already exists"



folder_name = folder
base_name = base_name
repeats = int(repeats)

original_fasta = fasta 

#for task in range(0, repeats+1) :
try :    
    tasks = [int(os.environ["SGE_TASK_ID"])-1]
    test = False
except KeyError :
    print "Not a cluster job?"
    tasks = range(0, repeats+1)
    test = False

for task in tasks:
    if task==0 :
        #Real
        input_file = original_fasta
        output_file = "%s/%s_palin_real.txt" \
                     % (folder_name, base_name)
    else :
        #Random
        input_file = "%s/random/random_%s_N%02i.fna" \
                     % (folder_name, base_name, task)
        output_file = "%s/random/%s_palin_random_N%02i.txt" \
                      % (folder_name, base_name, task)
        if not os.path.isfile(input_file) :
            print "Missing %s" % input_file
            write_shuffled_genome(original_fasta, input_file)
            print "Created it"

    print "Input: %s" % input_file
    print "Output: %s" % output_file
    start_time=time.time()
    genome_sequence = load_genome_sequence(input_file)
    #if test :
        #print ("this is in test mode and will not test the whole sequences")
        #genome_sequence = genome_sequence[:500] #For testing
    probability_threshold= 1.0 / len(genome_sequence)
    
    ################ actual running from here (min_pallindrome_len, max_pallindrome_len - options)
    palindrome_list = iteration_length(genome_sequence, min_pallindrome_len,max_pallindrome_len,\
                                       probability_threshold)
    write_palindromes_to_file(output_file, mask_NNN, palindrome_list)
    end_time=time.time()
    print 'Task %i took %.3f' %(task, end_time - start_time)
