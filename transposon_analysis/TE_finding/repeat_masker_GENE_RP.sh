cd /home/peter/Desktop/R.padi_Braker

#repeat modeler
#make the databsase first
# <RepeatModelerPath>/BuildDatabase -name Rpadi Rpadi_geneome.fa

/home/peter/Downloads/RepeatModeler/BuildDatabase -name Rpadi Rpadi.fa

RepeatModeler -pa 2 -database Rpadi 

make_repeatmodeler_database(name='a_database',input_filename='Rpadi.fa')

# -pa num_process -q is quick -s is slow -lcambig Outputs ambiguous DNA transposon fragments using a lower case name

#-source Includes for each annotation the HSP "evidence". -html 

# -u Creates an additional annotation file not processed by ProcessRepeats

RepeatMasker -pa 2 -u -qq -html -gff -source -lcambig -dir RP_repeat_masker_quick -species "Rhopalosiphum padi" ../R.padi_final_genome.v1_rename.fasta
 


