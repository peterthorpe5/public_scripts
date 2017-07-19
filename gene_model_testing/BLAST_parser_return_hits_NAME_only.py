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
    f_result = open(outfile + ".result.txt", 'w')
    f_result.write('#Query\tmismatches\tevalue\tpercent_identity\tdatabase\n')
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        alignment_hits = set([])
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    # For mismatches use (hsp.align_length - hsp.identities)
                    mmatches = hsp.align_length - hsp.identities
                    if mmatches <= mismatches:
                        f.write('Sequence: %s\n' % alignment.title)
                        f.write('Query_name = %s\n' %  blast_record.query)
                        f.write('query_length: %d\n' %  blast_record.query_length)
                        f.write('subject_length: %d\n' %  alignment.length)
                        f.write('e value: %s\n' %  hsp.expect)
                        # Calculate % of the matched region only:
                        percentage_hit = (hsp.identities)/ float(blast_record.query_length) * 100
                        f.write("percent_identiy: %1d\n" % percentage_hit)
                        f.write("%s\n" % hsp.query)
                        alignment_missing = hsp.match.replace(" ", "X")
                        db_entry = alignment.title
                        f.write("%s\n" % (alignment_missing.replace("|", " ")))
                        f.write("%s\n" % hsp.match[0:])
                        f.write("%s\n\n" % hsp.sbjct[0:])
                        result = "%s\t%d\t%s\t%.1f\t%s\n" % (blast_record.query,
                                                             mmatches,
                                                             hsp.expect,
                                                             percentage_hit,
                                                             alignment.title)
                        f_result.write(result)

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

