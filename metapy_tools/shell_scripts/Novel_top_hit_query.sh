
cd $HOME/scratch/tree_health/Illumina_results_2017Nov_first_run/version0.003

# this is a list of folder with control results in
folders=N01*_RESULTS
PROGAM=Swarm_d1

source ~/misc_python/THAPBI/THAPBI-pycits/venv-THAPBI-pycits-py3/bin/activate

export BLASTDB=/mnt/scratch/local/blast/ncbi/
# iterate through the results
for folder in ${folders}
do
	cd ./${folder}
	echo "NOW RUNNING:  ${folder}"
	Sample_name=${folder%????????}
	cd ./${Sample_name}_${PROGAM}
	cd ./clusters_results/
	cd ./novel*
	rm *.txt
	pwd
	python ~/public_scripts/metapy_tools/bin/Novel_top_hits.py
	cd ../../../../
done
