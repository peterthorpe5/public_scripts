#full info at end of file

"""Run RepeatModeler

RepeatModeler will produce consensus sequeces representing clusters of denovo repeat sequences, partialy classified by RepeatMasker """

#########################################################
from TE import *

#step 1

make_repeatmodeler_database(name='d_database',input_filename='all_mitochondria.fasta')

run_repeatmodeler('d_database') 
#use the RepeatModeler keyword to specify the path to your executable

print "IM breaking here... put the result into the Censor online and save in Censor_results"

#######################################################
#step2

censor_classifications = parse_online_censor('Censor_results.txt')

print_online_censor(censor_classifications,'censor_classifications_file')
                  

put_censor_classification_in_repeatmodeler_lib('mito.consensi.fa.classified.fasta',
                                                censor_classifications,
                                                'consensi.fa.censor')

#########################################################
#step3

"""Run RepeatMasker using the denovo lib

Some contig names are too long to be parsed in RepeatMasker.
However it is possible to replace the names with aliases and have
the translations in a file using the first function in this section.
It is important to remember to use the aliased genome assembly in the
other programs as well, so that redundancies can be resolved.


"""
###code_sequence_ids('genome_assembly',
##                  #'names_translations',
##                  #'coded_genome_assembly',
##                  #'prefix_for_aliases')

genome = 'Pea_aphid_genome_alt_names.fasta'

run_repeat_masker(genome, lib = 'consensi.fa.censor') 



#The default engine is ncbi and the default lib is eukaryota. RepeatModeler should make a folder in the CWD with the file 'consensi.fa.classified' in it. 
























#############################################################################################################
#############################################################################################################
#############################################################################################################
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
