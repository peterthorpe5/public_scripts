##############################################################################
#Title: Get nt_seq of interest, correspoding protein seq from the trinnotate
#   database file. Given a list of transcript of interest
###############################################################################

#author: Peter Thorpe September 2015. The James Hutton Insitute, Dundee, UK.

import os
import sys

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.2"
    sys.exit(0)


usage = """Use as follows:

$ python get_seq_cds_pep_corrspdng_transcript.py -t transcriptome_file -p proteins -c cds_file -d trinotate_database.txt -b top_hits.tab -n list _of_transcripts_of_interest-o output_prefix


#-d
trinotate databses in text format.
#-n
transcripts of interest, e.g. diff expresion list ....
#-b
top blast hit for the predicted peptides vs whatever ....(NR, swissprot....) formatted as ("\t")

Query tophit description

you can generate this by:
cat Mp_F.blastp.top_blast_hits.out | cut -f1,2,13 > ../Mp_F_blast_hits.tab


Give this script the required options and it will return the transcript, nt and amino acid sequence for a list of transcripts (-n) of interest

"""

#imports
from optparse import OptionParser
import datetime
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

###############################################################################
# function
###############################################################################

def wanted_genes(name_file):
    """function to retunr a list of wanted genes from file.
    """
    wanted = open(name_file, "r")
    names = wanted.readlines()
    name_data = [line.rstrip() for line in names
              if line.strip() != "" if not line.startswith("#")]
    wanted.close()
    name_set = set([])
    for i in name_data:
        name_set.add(i)
    return name_set

def tabular_file_to_dict(filename, key_column, value_columns):
    """Loads a tabular file into a Python dict.

    key_column = which column to use as the dict key.
    value_columns = list of columns to use as the dict values
    """
    answer = {}
    handle = open(filename, "r")
    for line in handle:
        if line.startswith("#"):
            continue
        if not line.strip():
            continue #if the last line is blank
        database_file_line = line.rstrip("\n").split("\t")
        value = ""
        for i in value_columns:
            value = value+database_file_line[i]
        key = database_file_line[key_column]
        answer[key] = value
    handle.close()
    return answer

def top_blast_hit_database(blast_output_file):
    return tabular_file_to_dict(blast_output_file, 0, [2])

#def tabular_file_to_info(trinotate):
    #return tabular_file_to_dict(trinotate, 1, [10, 11, 9])

def transcript_to_protein(trinotate):
    return tabular_file_to_dict(trinotate, 1, [5])


def tabular_file_to_info(trinnotate):
    handle = open(trinnotate, "r")
    transcript_info_dict = {}
    
    for line in handle:
        signal_p = ""
        TMMD = ""
        if line.startswith("#"):
            continue
        if not line.strip():
            continue #if the last line is blank
        database_file_line = line.rstrip("\n").split("\t")        
        transcript = database_file_line[1]
        # Record the desired description in the dict as a string,
        if database_file_line[10] !=".":
            signal_p = "Signal_p = Yes"
        if database_file_line[11] !=".":
            TMMD = "TMMD = Yes"
        if signal_p != "" and TMMD != "TMMD = Yes":
            info_out = "\tsecreted =  %s %s" %(signal_p, TMMD)
        else:
            info_out = " not_secreted"
            
        transcript_info_dict[transcript] = info_out
    return transcript_info_dict


def seq_getter(transcripome, proteins, cds, blast, trinotate, names_file, Output_prefix):
    """function to get the nt_seq/pep of genes of interest (names)- then opens up the
    annotation databse file and get the corresponsing protein seq which was predicted by
    transdecoder."""

    wanted = wanted_genes(names_file)
    ##########################################################################################                
    #   This is the trinnotate sql database in tab format
    #COLOUM 0 IS THE GENE/ COMPONENT ID
    #COLOUMN 1 IS THE TRANSCRIPT ID
    #COLOUMN 5 IS THE PROTEIN ID
    #gene_id	transcript_id	sprot_Top_BLASTX_hit	TrEMBL_Top_BLASTX_hit	RNAMMER	prot_id	prot_coords	sprot_Top_BLASTP_hit	TrEMBL_Top_BLASTP_hit	Pfam	SignalP	TmHMM	eggnog	gene_ontology_blast	gene_ontology_pfam	transcript	peptide
    ##########################################################################################
    
    transcript_info_dict = tabular_file_to_info(trinotate)
    transcript_to_protein_dict = transcript_to_protein(trinotate)
    top_blast_hit_dict = top_blast_hit_database(blast)
                             
    #nucleotide seq out
    nucleotide_out_file = Output_prefix+".nt.fasta"
    protein_out_file = Output_prefix+"_cds.pep.fasta"
    cds_nt_file_out = Output_prefix+"_cds.nt.fasta"
    #nt out
    f_nt_out = open(nucleotide_out_file, 'w')
    #protein file out
    f_PROTEIN_out = open(protein_out_file, 'w')
    #cds_file out
    f_cds_out = open(cds_nt_file_out, 'w')
    #index the fasta files
    transcriptome_record_db = SeqIO.index(transcripome, "fasta")
    cds_record_db = SeqIO.index(cds, "fasta")
    pep_record_db = SeqIO.index(proteins, "fasta")

    for i in wanted:
        if i in transcriptome_record_db:
            record = transcriptome_record_db[i]
            SeqIO.write(record, f_nt_out, "fasta")
            transdecoder_protein = transcript_to_protein_dict[i]
            
            if transdecoder_protein != ".":
                peprecord = pep_record_db[transdecoder_protein]
                try:
                    peprecord.description = top_blast_hit_dict[transdecoder_protein]+"\t"+transcript_info_dict[i]
                except KeyError:
                    # Join the three fields with spaces into one string:
                    peprecord.description = transcript_info_dict[i]#" ".join(transcript_info_dict[i])
                SeqIO.write(peprecord, f_PROTEIN_out, "fasta")
                
                cds_of_interest = cds_record_db[transdecoder_protein]
                cds_of_interest.description = peprecord.description

                # TODO - fill in descr
                SeqIO.write(cds_of_interest, f_cds_out, "fasta")

    f_nt_out.close()
    f_cds_out.close()
    f_PROTEIN_out.close()
    return True



##########################################################################################                
           
parser = OptionParser(usage=usage)

parser.add_option("-t", "--transcripome", dest="transcripome", default=None,
                  help="the transcriptome assembly",
                  metavar="FILE")

parser.add_option("-p", "--proteins", dest="proteins", default=None,
                  help="the predicted amino acid seq from the transcriptome assembly",
                  metavar="FILE")

parser.add_option("-c", "--cds", dest="cds", default=None,
                  help="the predicted cds from the transcripome assembly",
                  metavar="FILE")
parser.add_option("-d", "--trinotate", dest="trinotate", default=None,
                  help="the text file of the trinotate database,"
                  " used to get the corresponding cds/ pep for the transcript names",
                  metavar="FILE")

parser.add_option("-b", "--blast", dest="blast", default=None,
                  help="top hit descriptions for pep blasts",
                  metavar="FILE")

parser.add_option("-n", "--names", dest="names", default=None,
                  help="a text file list of transcript name of interest."
                  " This will be used to get the corressponding cds/pep",
                  metavar="FILE")

parser.add_option("-o", "--out", dest="Output_prefix", default=None,
                  help="Output_prefix filename",
                  metavar="FILE")

(options, args) = parser.parse_args()


def add_names_together(transcripome, filename):
    "function to get the default transdecoder names"
    if filename == "cds":
        return transcripome+".transdecoder.cds"
    if filename == "proteins":
        return transcripome+".transdecoder.pep"
    
#-t
transcripome = options.transcripome
#-p
proteins = options.proteins
#-c
cds = options.cds
#-b
blast = options.blast
#-d
trinotate = options.trinotate
#-n
names_file = options.names
# -o
Output_prefix = options.Output_prefix


#############

if not os.path.isfile(transcripome):
    sys_exit("Input transcriptome file not found: %s" % transcriptome)

if not os.path.isfile(proteins):
    sys_exit("Input proteins file not found: %s" % proteins)

seq_getter(transcripome, proteins, cds, blast, trinotate, names_file, Output_prefix)

   


print 'done'


