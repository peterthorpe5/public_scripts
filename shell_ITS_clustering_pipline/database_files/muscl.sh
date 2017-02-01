cd /home/pt40963/misc_python/THAPBI/Phyt_ITS_identifying_pipeline/database_files

 muscle -in Santi_database.trimed240.fasta -out Santi_database.trimed240.fasta.ali.fa
 
 muscle -in Santi_database.trimed240.fasta.ali.fa -out Santi_database.trimed240.fasta.refinefa -refine
 
 rm Santi_database.trimed240.fasta.ali.fa
 
 
 