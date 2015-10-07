"""Title: Pea_aphid_cds_finder

- basically the aphidbase cds- DNA file contains 3'
# and 5' regions. NOT the CDS of the predicted proteins.
#So, I generated a new mRNA file using GFFread. In the new
# seq_record.description() there are co-ordinates for the CDS. Are they correct?
This scripts pulls out these CDS and compares the tran;ated seq to the predicted amoni acid
seq. Only keeps it if they are the same.

"""

#biopython imports
from Bio.Seq import Seq

from Bio import SeqIO
import time

#os imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser


##################################################################################################################
##################################################################################################################
#functions

def CDS_translation(CDS):
    """find the CDS that starts with M by checking the 3 positive frames. If it cant find
    an M, then it defaults to fram +1"""
    protein = CDS.translate()
    if protein.startswith('M'):
        return protein, CDS
    nextframe = CDS[1:]
    protein2 = nextframe.translate()
    if protein2.startswith('M'):
        return protein2, CDS[1:]
    nextframe1 = CDS[2:]
    protein3 = nextframe1.translate()
    if protein3.startswith('M'):
        return protein3, CDS[2:]
    else:
        return protein, CDS # if the three positive frame have no meth start then retunr the first frame translation




def DNA_indexer(DNA_sequences):
    """ function to pull out the supposed CDS seq from the whole DNA_seq using
    the coordingates stated in the output for GFFread (using the GFF3 file
    from aphid base)"""
    # dict for the DNA specified by the CDS coordinates
    DNA_coordingate_dictionary= {}
    #put the new amino acid seq from the coordinate specified region
    DNA_coordinate_protein_seq_dic = {}
    
    for seq_record in SeqIO.parse(DNA_sequences, "fasta"):
        #some of the sequences have no coordinates for the CDS
        #only pick those that do
        if "CDS=" in seq_record.description:
            #split the description up until we get the info
            cds_cordinates = seq_record.description.split("CDS=")[1]
            start_coordinates = int(cds_cordinates.split("-")[0])
            end_coordinates = int(cds_cordinates.split("-")[1])
            new_CDS = seq_record.seq[(start_coordinates-3+1):end_coordinates]
        else:
            new_CDS = seq_record.seq
        #print new_CDS
        #call the function CDS_translation
        proteins_seq, CDS_starts_with_methionine = CDS_translation(new_CDS)
        #to keep the names the same in the DNA and predicted protein file
        #replace RA with PA - make life easier later
        new_name = seq_record.id.replace("RA", "PA")
        DNA_coordingate_dictionary [new_name] = CDS_starts_with_methionine
        DNA_coordinate_protein_seq_dic [new_name] = proteins_seq
    return DNA_coordingate_dictionary, DNA_coordinate_protein_seq_dic
        



def Pea_aphid_cds_finder(predicted_protein_sequences, DNA_sequences, outfile):
    """ this is the function that compares the translated NEW_CDS against the
    predicted proteins ... are they the same??? If so, write to file"""
    
    f_out= open(outfile, "w")
    #make a dictionary of the sequences. 
    #predicted_proteins =  SeqIO.index(predicted_protein_sequences, "fasta")
    #DNA_sequence_raw_with3_5UTR =  SeqIO.index(DNA_sequences, "fasta")
    
    #calls the function to pull out the CDS as defined by the coordinates
    DNA_coordingate_dictionary, DNA_coordinate_protein_seq_dic = DNA_indexer(DNA_sequences)
    yes_count = 0
    
    for seq_record in SeqIO.parse(predicted_protein_sequences, "fasta"):
        protein_seq = seq_record.seq
        gene_name = seq_record.id
        seq_record = DNA_coordinate_protein_seq_dic[seq_record.id]
        #print seq_record
        DNA_coordinate_protein = seq_record
        if str(protein_seq) == str(DNA_coordinate_protein):
            assert len(protein_seq) == len(DNA_coordinate_protein)
            DNA_coordinate_CDS = DNA_coordingate_dictionary[gene_name]
            # only keep the sequences which divide by three as the align back translate tool cries otherwise.
            if len(DNA_coordinate_CDS)%3 == 0:
                assert len(DNA_coordinate_CDS)%3 != 1
                assert len(DNA_coordinate_CDS)%3 != 2
                print >> f_out, ">%s\n%s" % (gene_name, DNA_coordinate_CDS)
                yes_count = yes_count+1
                #SeqIO.write(seq_record, f_out, "fasta")
    print yes_count
    return yes_count

######################################################################################################################################
####################################################################################################################################


#########################
#OK so the OptionParser doesnt work.... I find out how to call my function with it
#so this method until it is fixed




#Pea_aphid_cds_finder(argv[1],argv[2], argv[3])






if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

$ Pea_aphid_cds_finder -p predicted_protein_sequences -d DNA_sequences -o outfile



this will compare the new CDS/translated against the published pea aphid predicted proteins.
"""

parser = OptionParser(usage=usage)

parser.add_option("-p", dest="protein_file", default=None,
                  help="protein_file")
parser.add_option("-o", "--output", dest="out_file", default=None,
                  help="Output filename",
                  metavar="FILE")
parser.add_option("-d",  dest="DNA_file", default=None,
                  help="DNA_infile")



if len(args) > 1:
    stop_err("Expects no argument, one input filename")

(options, args) = parser.parse_args()

Pro_in_file = option.protein_file
out_file = options.out_file
DNA_file = options.DNA_file

(options, args) = parser.parse_args()

Pea_aphid_cds_finder(Pro_in_file, DNA_file, out_file)
