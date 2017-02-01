#!/usr/bin/env python
#author: Peter Thorpe September 2016. The James Hutton Insitute,Dundee,UK.

#Title:
#script to return ITS fasta regions from BLAST hits (tabular) "
# The BLAST should already have been perfomed:

#imports 
import os
from sys import stdin,argv


def split_blast_line(line):
    """function to split the blast line into its components"""
    queryId, subjectId, percIdentity, alnLength, mismatchCount, \
             gapOpenCount, queryStart, queryEnd, subjectStart, \
             subjectEnd, eVal, bitScore = line.split("\t")
    return queryId, subjectId, percIdentity, alnLength, mismatchCount, \
             gapOpenCount, queryStart, queryEnd, subjectStart, \
             subjectEnd, eVal, bitScore

def seq_getter(filename_in, blast, outfile):
    "fnction to get seq of interest"
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    f= open(outfile, 'w')
    blast_file = open(blast, "r")
    print ("indexing the genome")

    genome_databse = SeqIO.index(filename_in, "fasta")
    for line in blast_file:
        queryId, subjectId, percIdentity, alnLength, mismatchCount, \
             gapOpenCount, queryStart, queryEnd, subjectStart, \
             subjectEnd, eVal, bitScore = split_blast_line(line)
       
        seq_record =  genome_databse[subjectId.rstrip()]
        #print "yes"
        if int(subjectStart)> int(subjectEnd):
            # need to reverse complemnt
            seq_record.seq = seq_record.seq[int(subjectEnd): int(subjectStart)]
            seq_record.seq = seq_record.seq.reverse_complement()
            
        else:
            seq_record.seq = seq_record.seq[int(subjectStart): int(subjectEnd)]
        SeqIO.write(seq_record, f, "fasta")
        
    f.close()
    blast_file.close()
    return True


seq_getter(argv[1],argv[2], argv[3])


print 'done'

