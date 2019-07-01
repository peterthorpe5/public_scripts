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
   
upstream_of_TSS=" 20 200 500 1000 1500 2000 2500 3000 "   
downstream_of_TSS=" 0 20 200 500 1000 " 

for upstream in ${upstream_of_TSS}
do
    for downstream in ${downstream_of_TSS}
    do
        pycmd="python /storage/home/users/pjt6/public_scripts/TransStart/GFF_to_fasta.py \
        --gff MINC_TranStartbased_on_min_value.gff \
        -g Meloidogyne_incognita_V3_20131230_genome.fna \
        -x 4000 -i ${upstream} -u ${downstream} \
        -o MINC_TranStartbased_on_min_value_${downstream}_downstream_${upstream}_upstreamTSS.fasta
        --genome_gff sorted.final.gff"
        echo ${pycmd}
        eval ${pycmd}
    done
done
   
