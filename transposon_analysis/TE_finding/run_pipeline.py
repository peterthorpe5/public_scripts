#!/usr/bin/env python
# full info at end of file
from TE import *

"""Run RepeatModeler

RepeatModeler will produce consensus sequeces representing
clusters of denovo repeat sequences, partialy classified
by RepeatMasker
"""

#########################################################

# Make sure all prgrams are installed.
# You may need to rename the scaffold as this will break repeatmodeller
# Use the rewrite_as_fasta.py
print ("did you rename your scaffolds?")

#########################################################
# step 1

make_repeatmodeler_database(name='a_database',
                            input_filename='Rp.v1_alt.fasta')
print ("make_repeatmodeler_database done")

print ("...run_repeatmodeler")
run_repeatmodeler('a_database') 
# use the RepeatModeler keyword to specify the path to your executable

print ("""IM breaking here... put the result from within the
repeatmodeller directory: consensi.fa (or a_database)
into the Censor online  and save in Censor_results
 http://www.girinst.org/censor/
 """)

#########################################################
# step2

censor_classifications = parse_online_censor('Censor_results')

print_online_censor(censor_classifications,
                    'censor_classifications_file')

put_censor_classification_in_repeatmodeler_lib('consensi.fa.classified',
                                                censor_classifications,
                                                'consensi.fa.censor')
# step 2B
print ("""
 If you do this with multiple species, e.g. all the aphids, like I did.
 You will need to cat all the censor files and
 run cdhit at 100% to reduce the redundacy of all your
consensi.fa.censor files.""")

#########################################################
#step3

"""Run RepeatMasker using the denovo lib

Some contig names are too long to be parsed in RepeatMasker.
However it is possible to replace the names with aliases and
have the translations in a file using the first function
in this section. It is important to remember to use the
aliased genome assembly in the other programs as well,
so that redundancies can be resolved.

"""
code_sequence_ids('genome_assembly',
                  'names_translations',
                  'coded_genome_assembly',
                  'prefix_for_aliases')

genome = ""
#can just put genome seq in instead of coded genome seq. 

print ("... running run_repeat_masker")
run_repeat_masker('coded_genome_assembly',
                  lib = 'consensi.fa.censor') 



# The default engine is ncbi and the default lib is eukaryota.
# RepeatModeler should make a folder in the CWD with the file
# 'consensi.fa.classified' in it. 

#########################################################
# step 4
# once code section run - not via python:

 
/home/pt40963/Downloads/build_dictionary.pl --rm ./Pea_aphid_genome_alt_names.fasta.preThuMar171410582016.RMoutput/ --fuzzy > pea_fuzzy.txt
 
/home/pt40963/Downloads/OneCodeToFindThemAll.pl --dir /home/pt40963/Pea_aphid/altered_names/Pea_aphid_genome_alt_names.fasta.preThuMar171410582016.RMoutput --rm ./Pea_aphid_genome_alt_names.fasta.preThuMar171410582016.RMoutput/Pea_aphid_genome_alt_names.fasta.out --ltr pea_fuzzy.txt --fasta 



#step 5:
	   
	   # \# puting RepeatMasker results as parsed by OneCode... in the data structure
print ("parsing OCFA csv files")
TEs_from_OCFA, serial = parse_ocfa_elem_stored('/home/pt40963/Pea_aphid/altered_names/Pea_aphid_genome_alt_names.fasta.preThuMar171410582016.RMoutput/')
#
print ("DONE parsing OCFA csv files")

# \# adding LTRharvest results to the data structure
print ("reduce LTR into those")

TEs, serial = integrate_ltrharvest_to_RM_TEs('/home/pt40963/Pea_aphid/altered_names/pea_LTR_harvest.out',      
                                              '/home/pt40963/Pea_aphid/altered_names/Pea_aphid_genome_alt_names.fasta',  
                                              TEs_from_OCFA,
                                             serial)

print ("reduce PSI into those")
											  
#                                                                        
# \# adding TransposonPSI results to the data structure                                                                 
TEs = integrate_TransposonPSI_to_RM_TEs('/home/pt40963/Pea_aphid/altered_names/PSI_transposon_out/pea_alt.fasta.TPSI.allHits.chains.bestPerLocus', 
                                         '/home/pt40963/Pea_aphid/altered_names/Pea_aphid_genome_alt_names.fasta', 
                                         TEs, 
                                         serial)

print ("writing GFF")
# If this is set to true it break at score = max([int(s) for s in score.split('/')])    - cant int something in the line.
#unknow error for now.
write_gff(TEs, 'output.gff3', max_RM_OC_score=False)


############################################################################################################
# general info

"""README
TE discovery in a genome assembly
Overview

These are a set of wrappers compossing a pipline which finds transposable elements in a genome assembly. The pipline includes the following steps:

    MAKE A DENOVO LIB WITH REPEAT-MODELER. Input: genome assembly, Output: a library containing partially classified consensus sequences of de-novo repeat clusters. Program: RepeatModeler and all its many dependencies.

    ADD CLASSIFICATIONS FROM THE ONLINE CENSOR TEXT OUTPUT TO THE REPEAT LIB. Input: A) a library containing partially classified consensus sequences of de-novo repeat clusters. B) text file containing the online censor results. Output: a library containing more classified consensus sequences of de-novo repeat clusters (still some unknow classifications)

    SEARCH FOR TEs IN THE GENOME ASSEMBLY WITH REPEAT-MASKER USING THE DENOVO LIB. Input: genome assembly. Output: text file with a .out suffix containing classified TE loci with some redundancies. Program: RepeatMasker, which is a dependency of RepeatModeler.

    ELIMINATE REDUNDANCIES IN THE REPEAT-MASKER RESULTS. Input: text file with a .out suffix containing classified TE loci with some redundancies, or a directory with several .out files. Output: elem_stored.csv files, one per contig. Also generate other files per contig. Puts them in the same folder as the .out files used as input. Program: One Code To Find Them All.

    INDEPENDENT SEARCH FOR LTR ELEMENTS BASED ON SECONDARY STRUCTURE. Input: genome assembly, Output: text file with loci. Program: LTRharvest

    INDEPENDENT SEARCH FOR ELEMENTS BASED ON CODING SEQUENCES. Input: genome assembly, Output: text file with loci. Program: TransposonPSI

    ELIMINATE REDUNDANCIES AMONG PROGRAMS

        Read OneCodeTo... results. Input: Directory containing elem_stored.csv files. Output: Dictionary, Integer: num of elements in Dictionary.

        Read LTRharvest results and eliminate redundancies chosing the longer match between the programs. Input: Dictionary, Integer: num of elements in Dictionary. Output: Dictionary, Integer: num of elements in Dictionary.

        Read TransposonPSI results and eliminate redundancies chosing the longer match between the programs. Input: Dictionary, Integer: num of elements in Dictionary. Output: Dictionary, Integer: num of elements in Dictionary.

    PRINT NON-REDUNDANT OUTPUT FILE. Input: Dictionary. Output: gff3 file.


"""
