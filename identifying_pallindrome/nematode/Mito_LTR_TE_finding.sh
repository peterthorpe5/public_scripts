#$ -l hostname="n13*"
#$ -cwd

working_path=$HOME/nematode/Mitochondria_20160628/TE_finding

cd ${working_path}/

mkdir LTRharvest

cd LTRharvest

fasta_filenames=${working_path}/*.fa

# LTRharvest in genome tools  http://genometools.org/documents/ltrharvest.pdf
for f in ${fasta_filenames}
do
	echo "Running  ${f}"
	#step1
	#make the enhanced database
	cmd="gt suffixerator -db ${f} -indexname ${f} -tis -suf -lcp -des -ssp -sds -dna" 
	echo ${cmd}
	eval ${cmd}
	wait
	#step2	#search for the beasties. 
	cmd2="gt ltrharvest -index ${f} -mintsd 5 -maxtsd 100 > ${working_path}/TE_finding/LTRharvest/${f}_LTR_harvest.out" 
	echo ${cmd2}
	eval ${cmd2}
	
done

# (optional sequence clustering of output
#mkvtree .... look at the pdf. 

#############################################################################################
# use repeate modeller and repeat masker

cd ${working_path}/TE_finding

cd ${working_path}

cat *fa > ./TE_finding/all_mitochondria.fasta

cd ${working_path}/TE_finding
# this python script wraps the pipline. 
python run_pipeline.py