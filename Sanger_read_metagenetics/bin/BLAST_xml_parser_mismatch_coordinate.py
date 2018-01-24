#!/usr/bin/env python

# AUthor Peter Thorpe The James Hutton Institute Uk.
############################################################################
# Title: convert BLAST xml output to a custom format
############################################################################
"""
why? Just so it is easier to see the alignmnets.
"""
#############################################################################
# imports
import os
from sys import stdin, argv
import sys
from optparse import OptionParser
from Bio.Blast import NCBIXML
from collections import defaultdict

CONTROL_MIS_POSIT_COUNT = defaultdict(list)

def blast_file_opener(filename, evalue, mismatches, outfile):
    """Func takes in a BLAST xml output file (filename).
    writes out the various details of interests to the outfile.

    It filters the results based on evalue and number of
    mismatches, as defined by the user. """
    E_VALUE_THRESH = float(evalue)
    mismatches = int(mismatches)
    result_handle = open(filename)
    f = open(outfile, 'w')
    temp = outfile.split(".txt")[0]
    f_result = open(outfile.split(".txt")[0] + ".ALIGNMENT.txt", 'w')
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        alignment_hits = set([])
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    # For mismatches use (hsp.align_length - hsp.identities)
                    mmatches = hsp.align_length - hsp.identities
                    if mmatches <= mismatches:
                        mismatch_str = ""
                        alignment_missing = hsp.match.replace(" ", "X")
                        db_entry = blast_record.query
                        MisMatchX = 'X'
                        #print db_entry
                        mismatch_postions = [pos for pos, char in enumerate(alignment_missing) if char == MisMatchX]
                        for position in sorted(mismatch_postions):
                            #CONTROL_MIS_POSIT_COUNT[db_entry] = CONTROL_MIS_POSIT_COUNT[db_entry].append(position)
                            mismatch_str = mismatch_str + str(position) + "\t"
                        out_data = "%s\t%s\n" % (db_entry, mismatch_str)
                        f.write(out_data)
                        alignment_out = "%s\t%s\n" % (db_entry, alignment_missing)
                        f_result.write(alignment_out)
                        #print mismatch_postions

    f_result.close()
    f.close()
    result_handle.close()
    return alignment_hits


if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.2")
    sys.exit(0)


usage = """Use as follows:

converts

$ python BLAST... -i in.xml -m mismatches -e evalue -o out.txt

also prints out the seq and query of the matches with less than (mismatches) mismatches
"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="in_file", default=None,
                  help="in xml")
parser.add_option("-e", dest="evalue", default=0.00001,
                  help="evalue threshold for filtering")
parser.add_option("-m", dest="mismatches", default=2,
                  help="tell me the name of those whih has less than "
                  "this mismatches")
parser.add_option("-o", dest="out", default=None,
                  help="Output filename",
                  metavar="FILE")
(options, args) = parser.parse_args()
in_file = options.in_file
evalue = options.evalue
mismatches = options.mismatches
out = options.out

if __name__ == '__main__':
    blast_file_opener(in_file, evalue, mismatches, out)

