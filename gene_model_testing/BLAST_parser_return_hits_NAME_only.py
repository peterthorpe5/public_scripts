
import os
from sys import stdin,argv
import sys
from optparse import OptionParser

def blast_file_opener(filename, evalue, outfile):
    "this is a function to open up a .xml file blast results, the out put of\
is all the unique hit"
    evalue = float(evalue)
    result_handle = open(filename)
    output_filename = (outfile)
    f= open(output_filename, 'w')
    from Bio.Blast import NCBIXML
    blast_records = NCBIXML.parse(result_handle)
    #this function loads the blast lines into a list
    for blast_record in blast_records:
        #blast_records = list(blast_records)
        #print >>f,  blast_records
        E_VALUE_THRESH = evalue
        alignment_hits =set([])
    #for i in blast_records:    
        for alignment in blast_record.alignments:
            #print >>f,  blast_record
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    #print >>f,  '****Alignment****'
                    print >>f,  'Sequence:', alignment.title.split(" ")[1]#[19:36]
                    #print 'Sequence:', alignment.title.split(" ")[1]#[19:36]

                    
                    print >> f, 'Query_name = ', blast_record.query

                    #print 'Query_name = ', blast_record.query

                    print >> f,'query_length:', blast_record.query_length
                    print >> f,'subject_length:', alignment.length
                    print >> f,'e value:', hsp.expect
                    percentage_hit = (hsp.identities)/float(alignment.length)*100
                    print >> f,"percent_identiy: %1d" %(percentage_hit )
                    print >> f, hsp.query[0:]
                    print >> f, hsp.match[0:]
                    print >> f, hsp.sbjct[0:],'\n'
                    
    return alignment_hits


def write_out(filename, alignment_hits):
    output_filename = (filename)
    print >>f,  output_filename
    f= open(output_filename, 'w')
    #f.write("blast_output_condensed")
    #print >> f, alignment_hits
    f.close()
    return alignment_hits


if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

converts

$ python BLAST... -i in.xml -o out.txt
"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="in_file", default=None,
                  help="in xml")
parser.add_option("-e", dest="evalue", default=0.00001,
                  help="evalue threshold for filtering")
parser.add_option("-o", dest="out", default=None,
                  help="Output filename",
                  metavar="FILE")




(options, args) = parser.parse_args()

in_file = options.in_file
evalue = options.evalue 
out = options.out



blast_file_opener(in_file, evalue, out)

print 'done'


