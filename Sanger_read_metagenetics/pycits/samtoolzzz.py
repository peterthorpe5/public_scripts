#!/usr/bin/env python
# -*- coding: cp1252 -*-
#
# class to filter .sam output
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import subprocess
# from collections import Counter
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class Samtoolzzz_Error(Exception):
    """Exception raised when cd_hit fails"""
    def __init__(self, message):
        self.message = message

class Samtoolzzz(object):
    """class for filtering and manipulating sam output
    AIM:
    To get the number of reads that either have perfect mataches,
    or the number that hit with x mismatches
    """
    def __init__(self, infile):
        self.infile = infile


    def open_file(self):
        """fucntion to open the file
        return list split on \n"""
        with open(self.infile, 'r') as fh:
            data = fh.read().split("\n")
        return data


    def parse_samfile(self):
        """function class to parse the samfile and get the info
        required.
        
        QNAME: Query name of the read or the read pair
        FLAG: Bitwise flag (pairing, strand, mate strand, etc.)
        RNAME: Reference sequence name
        POS: 1-Based leftmost position of clipped alignment
        MAPQ: Mapping quality (Phred-scaled)
        CIGAR: Extended CIGAR string (operations: MIDNSHP)
        MRNM: Mate reference name (‘=’ if same as RNAME)
        MPOS: 1-based leftmost mate position
        ISIZE: Inferred insert size
        SEQQuery: Sequence on the same strand as the reference
        QUAL: Query quality (ASCII-33=Phred base quality
        TAG: XN:i:mismatches
        
        """
        # do I need this line?
        data = open_file(infile)
        # default dict to count the number of reads that hit
        # that db entery
        db_hits_counter = defaultdict(int)
        db_hits = defaultdict(list)
        for line in data:
            if line.startswith("@"):
                continue
            elements = line.split("\t")
            qname = elements[0]
            flag = elements[1]
            rname = elements[2]
            pos = elements[3]
            mapq = elements[4]
            cigar = elements[5]
            mrnm = elements[6]
            mpos = elements[7]
            isize = elements[8]
            seqquer = elements[9]
            qual = elements[10]
            tags = elements[11:]
            # populate the dict with a value when we have a read hitting
            # it
            db_hits_counter[rname] +=1
            db_hits[rname].append(qname + "\t" + seqquer)
        return db_hits_counter, db_hits, #qname, flag, rname, pos, \
               #mapq, cigar, mrnm, mpos, isize, seqquer, qual, tags
    
    def write_out_sequences(self, db_hits, seq_db):
        """function to obtain the sequences for the db hit and the reads
        so the user can visulase them as an alignmnet
        rname, seqquer are produced from the function about.
        seq_db need to be passed to it.
        """
        for keys, vals in db_hits:
            f_out = open(key + ".fasta", "w")
            # database entery - need to be passed the seq_db
            seq_record = seq_db[key]
            Seq_record.description = ""
            SeqIO.write(seq_record, f_out_file, "fasta")
            for entry in vals:
                name, seq = entry.split("\t")
                seq_record.id = name
                seq_record.seq = seq
                SeqIO.write(seq_record, f_out_file, "fasta")
            f_out.close()
            
            

    def write_out_read_to_db(db_hits):
        """function to simply write out the read and db hit."""
        for keys, vals in db_hits:
        

            
        
            
        
