
cd /home/pt40963/R.padi_genome/final_genome/altered_name

/home/pt40963/Downloads/build_dictionary.pl --rm /home/pt40963/R.padi_genome/final_genome/altered_name/Rp.v1_alt.fasta.preThuMar171408532016.RMoutput --fuzzy > rp_fuzzy.txt

/home/pt40963/Downloads/OneCodeToFindThemAll.pl --dir /home/pt40963/R.padi_genome/final_genome/altered_name/Rp.v1_alt.fasta.preThuMar171408532016.RMoutput --rm /home/pt40963/R.padi_genome/final_genome/altered_name/Rp.v1_alt.fasta.preThuMar171408532016.RMoutput/ --ltr rp_fuzzy.txt --strict --fasta


wait

python run_pipeline.py

