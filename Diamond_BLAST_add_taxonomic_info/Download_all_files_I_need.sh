#!/bin/bash
#$ -cwd
set -euo pipefail
# script to convert download NR

echo " NR database could be around 80GB. Diamond DB will be around 70 GB
other files will be 20GB... Do you have space here? 
KILL the program if not"


cd $PWD

echo "I DONT RECOMMEND THIS APPROACH. USE AT YOUR OWN RISK!!! This is not tested due to data size"

echo "im downloading files from NCBI"



wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz.md5
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
md5sum -c prot.accession2taxid.gz.md5
gunzip prot.accession2taxid.gz

wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxcat.zip
unzip taxcat.zip

wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz

wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.00.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.00.tar.gz.md5
md5sum -c nr.00.tar.gz.md5 < nr.00.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.01.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.01.tar.gz.md5
md5sum -c nr.01.tar.gz.md5 < nr.01.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.02.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.02.tar.gz.md5
md5sum -c nr.02.tar.gz.md5 < nr.02.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.03.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.03.tar.gz.md5
md5sum -c nr.03.tar.gz.md5 < nr.03.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.04.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.04.tar.gz.md5
md5sum -c nr.04.tar.gz.md5 < nr.04.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.05.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.05.tar.gz.md5
md5sum -c nr.05.tar.gz.md5 < nr.05.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.06.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.06.tar.gz.md5
md5sum -c nr.06.tar.gz.md5 < nr.06.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.07.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.07.tar.gz.md5
md5sum -c nr.07.tar.gz.md5 < nr.07.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.08.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.08.tar.gz.md5
md5sum -c nr.08.tar.gz.md5 < nr.08.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.09.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.09.tar.gz.md5
md5sum -c nr.09.tar.gz.md5 < nr.09.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.10.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.10.tar.gz.md5
md5sum -c nr.10.tar.gz.md5 < nr.10.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.11.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.11.tar.gz.md5
md5sum -c nr.11.tar.gz.md5 < nr.11.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.12.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.12.tar.gz.md5
md5sum -c nr.12.tar.gz.md5 < nr.12.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.13.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.13.tar.gz.md5
md5sum -c nr.13.tar.gz.md5 < nr.13.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.14.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.14.tar.gz.md5
md5sum -c nr.14.tar.gz.md5 < nr.14.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.15.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.15.tar.gz.md5
md5sum -c nr.15.tar.gz.md5 < nr.15.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.16.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.16.tar.gz.md5
md5sum -c nr.16.tar.gz.md5 < nr.16.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.17.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.17.tar.gz.md5
md5sum -c nr.17.tar.gz.md5 < nr.17.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.18.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.18.tar.gz.md5
md5sum -c nr.18.tar.gz.md5 < nr.18.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.19.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.19.tar.gz.md5
md5sum -c nr.19.tar.gz.md5 < nr.19.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.20.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.20.tar.gz.md5
md5sum -c nr.20.tar.gz.md5 < nr.20.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.21.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.21.tar.gz.md5
md5sum -c nr.21.tar.gz.md5 < nr.21.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.22.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.22.tar.gz.md5
md5sum -c nr.22.tar.gz.md5 < nr.22.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.23.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.23.tar.gz.md5
md5sum -c nr.23.tar.gz.md5 < nr.23.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.24.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.24.tar.gz.md5
md5sum -c nr.24.tar.gz.md5 < nr.24.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.25.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.25.tar.gz.md5
md5sum -c nr.25.tar.gz.md5 < nr.25.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.26.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.26.tar.gz.md5
md5sum -c nr.26.tar.gz.md5 < nr.26.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.27.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.27.tar.gz.md5
md5sum -c nr.27.tar.gz.md5 < nr.27.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.28.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.28.tar.gz.md5
md5sum -c nr.28.tar.gz.md5 < nr.28.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.29.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.29.tar.gz.md5
md5sum -c nr.29.tar.gz.md5 < nr.29.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.30.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.30.tar.gz.md5
md5sum -c nr.30.tar.gz.md5 < nr.30.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.31.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.31.tar.gz.md5
md5sum -c nr.31.tar.gz.md5 < nr.31.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.32.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.32.tar.gz.md5
md5sum -c nr.32.tar.gz.md5 < nr.32.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.33.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.33.tar.gz.md5
md5sum -c nr.33.tar.gz.md5 < nr.33.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.34.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.34.tar.gz.md5
md5sum -c nr.34.tar.gz.md5 < nr.34.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.35.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.35.tar.gz.md5
md5sum -c nr.35.tar.gz.md5 < nr.35.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.36.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.36.tar.gz.md5
md5sum -c nr.36.tar.gz.md5 < nr.36.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.37.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.37.tar.gz.md5
md5sum -c nr.37.tar.gz.md5 < nr.37.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.38.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.38.tar.gz.md5
md5sum -c nr.38.tar.gz.md5 < nr.38.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.39.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.39.tar.gz.md5
md5sum -c nr.39.tar.gz.md5 < nr.39.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.40.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.40.tar.gz.md5
md5sum -c nr.40.tar.gz.md5 < nr.40.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.41.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.41.tar.gz.md5
md5sum -c nr.41.tar.gz.md5 < nr.41.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.42.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.42.tar.gz.md5
md5sum -c nr.42.tar.gz.md5 < nr.42.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.43.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.43.tar.gz.md5
md5sum -c nr.43.tar.gz.md5 < nr.43.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.44.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.44.tar.gz.md5
md5sum -c nr.44.tar.gz.md5 < nr.44.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.45.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.45.tar.gz.md5
md5sum -c nr.45.tar.gz.md5 < nr.45.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.46.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.46.tar.gz.md5
md5sum -c nr.46.tar.gz.md5 < nr.46.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.47.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.47.tar.gz.md5
md5sum -c nr.47.tar.gz.md5 < nr.47.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.48.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.48.tar.gz.md5
md5sum -c nr.48.tar.gz.md5 < nr.48.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.49.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.49.tar.gz.md5
md5sum -c nr.49.tar.gz.md5 < nr.49.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.50.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.50.tar.gz.md5
md5sum -c nr.50.tar.gz.md5 < nr.50.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.51.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.51.tar.gz.md5
md5sum -c nr.51.tar.gz.md5 < nr.51.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.52.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.52.tar.gz.md5
md5sum -c nr.52.tar.gz.md5 < nr.52.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.53.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.53.tar.gz.md5
md5sum -c nr.53.tar.gz.md5 < nr.53.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.54.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.54.tar.gz.md5
md5sum -c nr.54.tar.gz.md5 < nr.54.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.55.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.55.tar.gz.md5
md5sum -c nr.55.tar.gz.md5 < nr.55.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.56.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.56.tar.gz.md5
md5sum -c nr.56.tar.gz.md5 < nr.56.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.57.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.57.tar.gz.md5
md5sum -c nr.57.tar.gz.md5 < nr.57.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.58.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.58.tar.gz.md5
md5sum -c nr.58.tar.gz.md5 < nr.58.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.59.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.59.tar.gz.md5
md5sum -c nr.59.tar.gz.md5 < nr.59.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.60.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.60.tar.gz.md5
md5sum -c nr.60.tar.gz.md5 < nr.60.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.61.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.61.tar.gz.md5
md5sum -c nr.61.tar.gz.md5 < nr.61.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.62.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.62.tar.gz.md5
md5sum -c nr.62.tar.gz.md5 < nr.62.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.63.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.63.tar.gz.md5
md5sum -c nr.63.tar.gz.md5 < nr.63.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.64.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.64.tar.gz.md5
md5sum -c nr.64.tar.gz.md5 < nr.64.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.65.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.65.tar.gz.md5
md5sum -c nr.65.tar.gz.md5 < nr.65.tar.gz



wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.66.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.66.tar.gz.md5
md5sum -c nr.66.tar.gz.md5 < nr.66.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.67.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.67.tar.gz.md5
md5sum -c nr.67.tar.gz.md5 < nr.67.tar.gz

wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.68.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.68.tar.gz.md5
md5sum -c nr.68.tar.gz.md5 < nr.68.tar.gz


wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.69.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nr.69.tar.gz.md5
md5sum -c nr.69.tar.gz.md5 < nr.69.tar.gz

echo "please check there should only be 69 of the nr files as of Jan 2023. This changes!!"


folder=*tar.gz

for fol in ${folder}
do
	tar -zxvf ${fol}
done

echo downloading and unzipping done

export BLASTDB=$`pwd` 

blastdbcmd -entry 'all' -db nr > nr.faa

echo "downloading and unzipping done"

python prepare_accession_to_description_db.py

