#$ -l hostname="n13*"
#$ -cwd

cd /home/pt40963/Pea_aphid/altered_names

# LTRharvest in genome tools  http://genometools.org/documents/ltrharvest.pdf

#make the enhanced database
gt suffixerator -db Pea_aphid_genome_alt_names.fasta -indexname pea -tis -suf -lcp -des -ssp -sds -dna


#search for the beasties. 

gt ltrharvest -index pea -mintsd 5 -maxtsd 100 > pea_LTR_harvest.out


# (optional sequence clusterin gof output

#mkvtree .... look at the pdf. 
