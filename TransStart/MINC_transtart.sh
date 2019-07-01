cd /shelf/apps/pjt6/mincognita

module load samtools

conda activate python36


#samtools sort -@ 2 -o minc_sorted.bam  minc.bam 



# transend
python /storage/home/users/pjt6/public_scripts/TransStart/TrannEnd_no_exon.py \
-b minc_sorted.bam \
-g Meloidogyne_incognita_V3_20131230_genome.fna \
--gff sorted.final.gff \
-o MINC_End.txt \
--keep_gene_depth yes

 python /storage/home/users/pjt6/public_scripts/TransStart/GFF_to_fasta.py \
 --gff MINC_Endbased_on_min_value.gff \
 -g Meloidogyne_incognita_V3_20131230_genome.fna \
 
 # transtart
python ./TransStart_no_exon.py \
-b minc_sorted.bam \
-g Meloidogyne_incognita_V3_20131230_genome.fna \
--gff sorted.final.gff \
-o MINC_TranStart.txt \
--keep_gene_depth yes

python /storage/home/users/pjt6/public_scripts/TransStart/GFF_to_fasta.py \
 --gff MINC_TranStartbased_on_min_value.gff \
 -g Meloidogyne_incognita_V3_20131230_genome.fna \
 
 
python /storage/home/users/pjt6/public_scripts/TransStart/GFF_to_fasta.py \
   --gff TransStart_walk3based_on_min_value.gff -g nGr.v1.1.fa \
   -o GROS_TSS_predition_maxLen4000_upstream1000_200_intoUTR.fasta \
   -x 4000 -i 200 -u 1000