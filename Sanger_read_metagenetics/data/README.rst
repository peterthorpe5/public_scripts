READme: generating the db from NCBI
==================================

Representative sequences for the beastie of interest was BLASTn searched against GenBank nt database,
evalue 0.0001 (relaxed).

Make a fasta of nt for searching through:
========================================
export BLASTDB=/mnt/scratch/local/blast/ncbi/
blastdbcmd -entry 'all' -db nt > nt.faa

 
Blast the representive against nt
=================================
export BLASTDB=/mnt/scratch/local/blast/ncbi/
blastn -db nt -query pythium_its.fasta -evalue 0.001 -dust no -num_threads 6 -outfmt 5 -out ITS1_vs_NR.xml

Parse the xml for info of interest
==================================
python ~/public_scripts/Sanger_read_metagenetics/bin/BLAST_xml_parser.py -i ITS1_vs_NR.xml -e 0.001 -m 20 -o Pythium ITS1_vs_NR.txt   

# may also need to cut to specific names for this
cut -f5 Pythium.result.txt > to_get.txt

get the sequences from the nt database generated above
======================================================
python ~/misc_python/get_sequences_i_want_from_fasta_command_line.py ~/scratch/blast_databases/nt.faa to_get.txt Pythium_genebank_18s_ITS1.fasta

