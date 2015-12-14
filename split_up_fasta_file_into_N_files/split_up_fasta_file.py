#!/usr/bin/env python

# title: script to split up a massive fasta file into smaller files for
#blasting multiple files rather than one large file.

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
from sys import stdin,argv
import sys
from optparse import OptionParser

############ imports

# make a folder for the outfiles


def make_me_a_folder():
    directory = os.getcwd()
    try:
        os.mkdir('split_fasta_files')
    except OSError:
        print "already exists"
    #os.chdir("split_fasta_files")
    return True

########################################
    
def count_size_of_fasta(in_fasta):
    "function to find out how many sequences"
    with open(in_fasta, "r") as fasta_file:
        count = 0
        for line in fasta_file:
            if line.startswith(">"):
                count = count+1
        return count

def how_many_seq_per_file(how_many_files, in_fasta):
    "function to find out how many seq in each file"
    total_size = count_size_of_fasta(in_fasta)
    how_many_seq_per_file = total_size/int(how_many_files)
    print "number of seq per file: ", how_many_seq_per_file
    return how_many_seq_per_file, total_size

def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.
 
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
 
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch
    
 
def seq_getter(how_many_files, filename):
    "function to write the outfiles"
    names_printed = set([])
    how_many_seq_in_file, total_size = how_many_seq_per_file(how_many_files, filename)
    record_iter = SeqIO.parse(filename, "fasta")
    
    for i, batch in enumerate(batch_iterator(record_iter, how_many_seq_in_file)) :
        outfile = "./split_fasta_files/"+filename.split("fa")[0]+"_%d.fasta" %(i+1)
        handle = open(outfile, "w")
        count = SeqIO.write(batch, handle, "fasta")
        handle.close()
        #print "Wrote %i records to %s" % (count, filename)

    return True


#################################################################################################

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

$ split_up_fasta_file.py -n how_many_files -i in.fasta

this is a script to splits a big fasta file up into smaller ones.

this program will make a folder called split_fasta_files and
put the fasta files in there

"""

parser = OptionParser(usage=usage)

parser.add_option("-n", "--num", dest="How_many_file", default=None,
                  help="How_many_file to split the fasta up into")

parser.add_option("-i", "--in", dest="in_file", default=None,
                  help="in_file.fasta",
                  metavar="FILE")




(options, args) = parser.parse_args()

How_many_file = options.How_many_file
in_file = options.in_file
make_me_a_folder


# run the fuctions
make_me_a_folder()
seq_getter(How_many_file, in_file)
