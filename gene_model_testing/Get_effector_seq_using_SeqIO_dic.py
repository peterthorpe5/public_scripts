##############################################################################
#Title: Get nt_seq of interest, correspoding protein seq effector balst matches
###############################################################################

#imports
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import os
from sys import stdin,argv
import sys
from optparse import OptionParser

###############################################################################

# function

###############################################################################

file_name = 'test.txt'
script_dir = os.getcwd()
dest_dir = os.path.join(script_dir, 'Bos2010')
try:
    os.makedirs(dest_dir)
except OSError:
    print "already exists"

######################################################################################################
def seq_getter(blast_hits_wanted, cds_file, protein_file):
    """this is a function to get the nt_seq of genes of interest and get the corresponsing
    protein seq which was predicted by transdecoder.

    ."""
    blast_hits_wanted = open(blast_hits_wanted, "r")
 
    ##########################################################################################                
    #effector list
    ##########################################################################################
    
    names_of_effectors = """Mp1_MpSec22
Mp2
Mp3
Mp4_MpSec74
Mp5
Mp6
Mp7
Mp8
Mp10_MpSec34
Mp11
Mp12
Mp14
Mp15_MpSec24
Mp16_MpSec39
Mp17
Mp19
Mp20
Mp21
Mp22
Mp23
Mp24
Mp28
Mp29
Mp30
Mp31
Mp32
Mp33
Mp35
Mp36
Mp37
Mp38
Mp39
Mp40
Mp41
Mp42_MpSec35
Mp43
Mp44
Mp45_MpSec50
Mp46
Mp47_MpSec26
Mp49
Mp50
Mp51
Mp53
Mp54
MpCOO2
MpSec1_noSP
MpSec2_NoSP_partial
MpSec7_partial
MpSec5
MpSec8
MpSec3_partial
MpSec9
MpSec25
MpSec26
MpSec27_FL
MpSec28
MpSec38
MpSec42
MpSec75
Dystrophin_DMEL""".split()
##    nhandles = dict()
##    for gene in names_of_effectors:
##        filename = "./Bos2010/"+gene+"_cds.fasta"
##        nhandles[gene] = open (filename, "w")
    phandles = dict()
    for gene in names_of_effectors:
        filename = "./Bos2010/"+gene+"_pep.fasta"
        phandles[gene] = open (filename, "w")
    
    ##########################################################################################                
    #actula code block
    ##########################################################################################
    print("Indexing...")
    names_already_printed = set([])
    Bos_2010_pep_seq = SeqIO.index("/home/pt40963/gene_model_testing/Mp_candidates_Bos_lab_march2014.fasta", "fasta")
    #Bos_2010_nt_seq = SeqIO.index("/home/pt40963/gene_model_testing/Bo2010_nt_proper_length.fasta", "fasta")

    protein_sequences =  SeqIO.index(protein_file, "fasta")
    #nucleotide_sequences =  SeqIO.index(cds_file, "fasta")
    print("Starting output...")
    for line in blast_hits_wanted:
        if line.startswith("#"):
            continue
        data = line.rstrip("\n").split("\t")
        gene = data[0]
        blast_hit_matches = data[1]

        #seq_record = nucleotide_sequences[blast_hit_matches]
        #SeqIO.write(seq_record, nhandles[gene], "fasta")

        seq_record = protein_sequences[blast_hit_matches]
        SeqIO.write(seq_record, phandles[gene], "fasta")
        if not gene in names_already_printed:        
            seq_record = Bos_2010_pep_seq[gene]
            SeqIO.write(seq_record, phandles[gene], "fasta")
            names_already_printed.add(gene)
        

##
##    for gene in names_of_effectors:
##        nhandles[gene].close()
    for gene in names_of_effectors:
        phandles[gene].close()
    return True

#################################################################################################

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

converts

$ python BLAST... -i in.xml -o out.txt
"""

parser = OptionParser(usage=usage)

parser.add_option("-b", dest="blast_output", default=None,
                  help="in xml")

parser.add_option("-p", dest="protein_file", default=None,
                  help="predicted amino acid seq of genes")
parser.add_option("-n", dest="nuc", default=None,
                  help="predcited nt cds of genes. ",
                  metavar="FILE")




(options, args) = parser.parse_args()

blast_output = options.blast_output
protein_file = options.protein_file
nuc = options.nuc 

seq_getter(blast_output, nuc, protein_file)


print 'done'

##################################################################################################



#seq_getter('all_up_bodies_hits.out', \
           #'all_aphids_cds.fasta',\
           #'all_aphids_proteins.fasta')
