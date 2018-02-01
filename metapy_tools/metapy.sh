cd $HOME//scratch/tree_health/Illumina_results_2017Nov_first_run/version0.003

#cp /mnt/shared/projects/Phytophthora/201711_ITS1_metabarcoding_1_THAPBI/static/RawData/* ./

source ~/misc_python/THAPBI/THAPBI-pycits/venv-THAPBI-pycits-py3/bin/activate

# run metapy. On a 96 samples this take about 2 day (without BLASTCLUST)
python generate_metapy.py
##########################################################################################
wait
python interogate_results.py Control_C1
python interogate_results.py Control_C2
python interogate_results.py Control_C3
python interogate_results.py Control_C4
python interogate_results.py Control

echo "python interogate_results.py austrocedri"
python interogate_results.py austrocedri
echo "python interogate_results.py ramorum"
python interogate_results.py ramorum
echo "python interogate_results.py lateralis"
python interogate_results.py lateralis
echo "python interogate_results.py cactorum"
python interogate_results.py cactorum
echo "python interogate_results.py gonapodyides"
python interogate_results.py gonapodyides
echo "python interogate_results.py agathidicida"
python interogate_results.py agathidicida


##########################################################################################

#mkdir reads

mv *fastq* ./reads

#  dont delete for now
#rm -rf reads

mkdir swarm_results
mkdir bowtie_results
mkdir $HOME/scratch/tree_health/Illumina_results_2017Nov_first_run/version0.003/Novel_blast_results_swarm


############################################################################################
folders=*_RESULTS
PROGAM=Swarm_d1

export BLASTDB=/mnt/scratch/local/blast/ncbi/
# iterate through the results
for folder in ${folders}
do
	cd ./${folder}
	echo "NOW RUNNING:  ${folder}"
	Sample_name=${folder%????????}
	cd ./${Sample_name}_${PROGAM}
	cp .RESULTS $HOME/scratch/tree_health/Illumina_results_2017Nov_first_run/version0.003/swarm_results
	cd ./clusters_results/ 
	cd ./novel*
	rm *.txt
	pwd
	python ~/public_scripts/metapy_tools/bin/Novel_top_hits.py
	cp *.txt $HOME/scratch/tree_health/Illumina_results_2017Nov_first_run/version0.003/Novel_blast_results_swarm
	cd ../../../../
done


############################################################################################
folders=*_RESULTS
PROGAM=bowtie


export BLASTDB=/mnt/scratch/local/blast/ncbi/
# iterate through the results
for folder in ${folders}
do
	cd ./${folder}
	echo "NOW RUNNING:  ${folder}"
	Sample_name=${folder%????????}
	cd ./${Sample_name}_${PROGAM}
	cp .RESULTS $HOME/scratch/tree_health/Illumina_results_2017Nov_first_run/version0.003/bowtie_results
	cd ./clusters_results/ 
	cd ./novel*
	cd ../../../../
done

########################################
# summerize:
cd bowtie_results
grep -H "percent of assembled-reads clustering with database" * > bowtie_percentage_cluster_with_db.txt
grep -H "Total number of assembled sequences" * > total_number_of_assembled_sequences.txt

cd ../
cd swarm_results
grep -H "percent of assembled-reads clustering with database" * > swarm_percentage_cluster_with_db.txt
grep -H "Total number of assembled sequences" * > total_number_of_assembled_sequences.txt


