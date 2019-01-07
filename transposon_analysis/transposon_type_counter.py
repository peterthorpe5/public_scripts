# Title: script to count transposon types
# author: Peter Thorpe September 2015. The James Hutton Insitute, Dundee, UK.

# imports
import random
import os
import sys
from optparse import OptionParser
import datetime
from collections import defaultdict


def parseTransposon_info(transposon_info):
    """function to parse the element of trasposon_file
    . The GFF output from the transposon finding pipeline
    is used here to define the coordinate and annot of the
    transposon. NOTE. The scaffold names may have been changed
    so they fit into the program prior to this. If this brakes.
    check them
    INPUT e.g.
    Mc2825	Repeatmasker-OneCode	transposable_element	26019	26484	2108	-	.	ID=element156452;Name=DNA/hAT-Tip100;Note=?
    Mc2825	Repeatmasker-OneCode	transposable_element	27906	27982	300	-	.	ID=element156453;Name=DNA/RC;Note=?
    Mc2825	Repeatmasker-OneCode	transposable_element	6290	6427	488	-	.	ID=element156450;Name=DNA/hAT;Note=?
    """
    note_annot = transposon_info.split(";Note=")[1]
    # eg.g ID=element96159;Name=DNA/MULE-MuDR;Note=?
    Name_annot = transposon_info.split("Name=")[1]
    Name_annot = Name_annot.split(";")[0]
    return Name_annot, note_annot.rstrip()


def parse_gff(gff):
    """function to parse GFF and produce a dictionary of transposon to count
    s...
    
    function to parse the element of trasposon_file
    . The GFF output from the transposon finding pipeline
    is used here to define the coordinate and annot of the
    transposon. NOTE. The scaffold names may have been changed
    so they fit into the program prior to this. If this brakes.
    check them
    INPUT e.g.
    Mc2825	Repeatmasker-OneCode	transposable_element	26019	26484	2108	-	.	ID=element156452;Name=DNA/hAT-Tip100;Note=?
    Mc2825	Repeatmasker-OneCode	transposable_element	27906	27982	300	-	.	ID=element156453;Name=DNA/RC;Note=?
    Mc2825	Repeatmasker-OneCode	transposable_element	6290	6427	488	-	.	ID=element156450;Name=DNA/hAT;Note=?
    """
    f_in = open(gff, "r")
    # data to output of function
    scaffold_to_gene_dict = dict()
    gene_to_exon_count = dict()
    gene_start_stop_dict = dict()
    count = 1
    transposon_count_dict = defaultdict(int)
    note_dict = defaultdict(int)
    # sum up the commumalitve lenght of a type of transpson in the genome
    # why? to see if certaqin TE take up way more space than they should
    # but that need more downstream analysis of course.
    TE_comm_len = defaultdict(int) # default value of int is 0
    TE_comm_len_name_annot = defaultdict(int) # default value of int is 0

    # iterate through gff
    for line in f_in:
        if line.startswith("#"):
            continue
        if not line.strip():
                continue  #  if the last line is blank
        scaffold, aug, cds_type, start, stop, e, f, \
                  g,gene_info = line.split("\t")
        Name_annot, note_annot = parseTransposon_info(gene_info)
        transposon_count_dict[Name_annot] += 1
        note_dict[note_annot] += 1
        TE_comm_len[Name_annot] += (int(stop) - int(start))
        # note annot is a much more refined type of TE 
        TE_comm_len_name_annot[note_annot] += (int(stop) - int(start))
    f_in.close()
    #print(TE_comm_len)
    #print(TE_comm_len_name_annot)
    return transposon_count_dict, note_dict, TE_comm_len, TE_comm_len_name_annot


if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.3")
    sys.exit(0)


usage = """Use as follows:

$ python transposon_type_counter.py -t transposon.gff -o outfile.out

script reports the types and number of transpsons
"""

parser = OptionParser(usage=usage)

parser.add_option("-t", "--gff", "--transposons",
                  dest="transposons_gff",
                  default=None,
                  help="gff file generated for the trapsosons finding",
                  metavar="FILE")
parser.add_option("-o", "--out",
                  dest="out_file",
                  default="transposon.counts",
                  help="outfile for the results",
                  metavar="FILE")


(options, args) = parser.parse_args()

#-o
out_file= options.out_file
#-t
transposons_gff= options.transposons_gff

#########################################################
# run the program
# Run as script
if __name__ == '__main__':
    # call the main function
    if not os.path.isfile(transposons_gff):
        print("sorry, couldn't open the file: " % transposons_gff)
    transposon_count_dict, note_dict, TE_comm_len, \
                           TE_comm_len_name_annot = parse_gff(transposons_gff)
    outfile = open(out_file, "w")

    title = "#Transposon\tCount\n"
    outfile.write(title)
    for key, vals in transposon_count_dict.items():
        outfmt = "%s\t%d\n" % (key, vals)
        outfile.write(outfmt)
    outfile.close()
    
    ###################################
    # NOTE is the "more refined" type of TE
    note_out = out_file.replace(".txt", "") + "_NOTES.txt"
    noteout = open(note_out, "w")
    title = "#Notes\tCount\n"
    noteout.write(title)
    for key, vals in note_dict.items():
        outfmt = "%s\t%d\n" % (key, vals)
        noteout.write(outfmt)
    noteout.close()

    ###################################
    # output the commumalative len of the TE in the genome
    note_out = out_file.replace(".txt", "") + "_commu_length_TE.txt"
    noteout = open(note_out, "w")
    title = "#Notes\ttotal_length\n"
    noteout.write(title)
    for key, vals in sorted(TE_comm_len.items()):
        outfmt = "%s\t%d\n" % (key, vals)
        noteout.write(outfmt)
    noteout.close()
    
    ###################################
    # output the commumalative len of the TE in the genome
    note_out = out_file.replace(".txt", "") + "_commu_length_Family.txt"
    noteout = open(note_out, "w")
    title = "#Notes\ttotal_length\n"
    noteout.write(title)
    for key, vals in sorted(TE_comm_len_name_annot.items()):
        outfmt = "%s\t%d\n" % (key, vals)
        noteout.write(outfmt)
    noteout.close()



