#!/bin/bash
#$ -cwd
 ### not needed   #$ -l hostname="n17*"
cd $HOME/scratch/tree_health/soils_scot_gov

# FILL THESE IN:
interest_list="Control_C1 Control_C2 Control_C3 Control_C4 Control austroc ramorum 
lateralis cactorum gonapodyides agathidicida"

Working_directory=$HOME/scratch/tree_health/soils_scot_gov
project="Soils_scot_gov"

Raw_illumna_data=/mnt/shared/projects/Phytophthora/201808_Soils_Scot_Gov_ITS_Metabarcoding/RawData/

# Nothing to fill in from here!
#######################################


# copy reads from source. Never alter the master copy.
 cp ${Raw_illumna_data}/* ./
rm *_I*_001.fastq.gz
rm Undetermined_*
# activate the VM
source ~/misc_python/THAPBI/THAPBI-pycits/venv-THAPBI-pycits-py3/bin/activate

# run metapy. On a 96 samples this take about 1 day (without BLASTCLUST and graph drawing)
# This py script calls metapy with all the read in a folder. 
python ~/public_scripts/metapy_tools/generate_metapy.py

wait
mkdir reads
mkdir swarm_results
mkdir bowtie_results
mkdir Novel_blast_results_swarm

mv *fastq* ./reads
#  dont delete for now
rm -rf reads


############################################################################################
# This nest of code blast the novel clusters against nt. Only those clusters which are over
# threshold in size. see ~/public_scripts/metapy_tools/bin/Novel_top_hits.py for more details
# Also copies the results to folders or the clustering and the BLAST
folders=*_RESULTS
PROGAM=Swarm_d1

export BLASTDB=/mnt/apps/databases/blast-ncbi/
# iterate through the results
for folder in ${folders}
do
	cd ./${folder}
	echo "NOW RUNNING:  ${folder}"
	Sample_name=${folder%????????}
	cd ./${Sample_name}_${PROGAM}
	echo 
	cp *.RESULTS ${Working_directory}/swarm_results
	cd ./clusters_results/ 
	cd ./novel*
	#rm *.txt
	pwd
	python ~/public_scripts/metapy_tools/bin/Novel_top_hits.py --min_cluster_size 20 --threads 16
	cp *.txt ${Working_directory}/Novel_blast_results_swarm
	cd ${Working_directory}
done

###############################################################################
# summerise species found:

cd ./swarm_results
for species in ${interest_list}
do
	echo "NOW getting:  ${species}"
	echo "#Filename	Species	Number_of_reads" > ${species}_swarm_identified.out
	grep -H "${species}" *.RESULTS >> ${species}_swarm_identified.out
done

cd ../

############################################################################################
# Also copies the results to folders of the clustering
folders=*_RESULTS
PROGAM=bowtie

export BLASTDB=/mnt/apps/databases/blast-ncbi/
# iterate through the results
for folder in ${folders}
do
	cd ./${folder}
	echo "NOW RUNNING:  ${folder}"
	Sample_name=${folder%????????}
	cd ./${Sample_name}_${PROGAM}
	cp *RESULTS ${Working_directory}/bowtie_results
	cd ./clusters_results/ 
	cd ./novel*
	cd ${Working_directory}
done

########################################
# summerize the % of reads clustered with the database:
cd bowtie_results
grep -H "percent of assembled-reads clustering with database" *.RESULTS > bowtie_percentage_cluster_with_db.txt
grep -H "Total number of assembled sequences" *.RESULTS > total_number_of_assembled_sequences.txt

cd ../
cd swarm_results
grep -H "percent of assembled-reads clustering with database" *.RESULTS > swarm_percentage_cluster_with_db.txt
grep -H "Total number of assembled sequences" *.RESULTS > total_number_of_assembled_sequences.txt

cd ${Working_directory}

######################################
# search the bowtie reasults for species of interest
cd ./bowtie_results
for species in ${interest_list}
do
	echo "NOW getting:  ${species}"
	echo "#Filename	Species	Number_of_reads" > ${species}_bowtie_identified.out
	grep -H "${species}" *.RESULTS >> ${species}_bowtie_identified.out
done

cd ${Working_directory}


cd ${Working_directory}/swarm_results
python ~/public_scripts/metapy_tools/results_summary.py -t 50  -o ${project}_Swarm_RESULTS

cd ${Working_directory}/bowtie_results/
python ~/public_scripts/metapy_tools/results_summary.py  -t 10  -o ${project}_Bowtie_RESULTS

cd ${Working_directory}/Novel_blast_results_swarm/
python ~/public_scripts/metapy_tools/results_summary.py --blast yes  -o ${project}_Novel_BLAST_results

echo "FINISHED"

