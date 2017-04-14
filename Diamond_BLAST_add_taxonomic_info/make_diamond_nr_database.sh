#!/bin/bash


# script to convert NR database to fasta files and make a diamond blastdatabse.


cd $PWD
# Fill these out


#################################################################
# Help with diamond
# to install diamond from source
#wget http://github.com/bbuchfink/diamond/archive/v0.7.9.tar.gz
#tar xzf v0.7.9.tar.gz
#cd diamond-0.7.9/src
# # optional, for installing Boost
##./install-boost 
#make
# 
 
#put the /bin in your PATH

cd /path_to/Downloads/blast_databases

export BLASTDB=/mnt/shared-new/cluster/blast/ncbi/extracted/

# can only use protein databases with this program.
blastdbcmd -entry 'all' -db nr > nr.faa

echo im making the nr fasta file

/path_to/Downloads/diamond-0.7.9/bin/diamond makedb --in nr.faa -d nr

echo nr fasta done
#diamond makedb --in /mnt/shared/cluster/blast/ncbi/extracted/uniprot_sprot.faa -d uniprot

#diamond makedb --in /mnt/shared/cluster/blast/ncbi/extracteduniref90.faa -d uniref90


#files required for pyhon script to get tax id and species name ..

echo im downloading files from NCBI

wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz.md5
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
md5sum -c prot.accession2taxid.gz.md5
gunzip prot.accession2taxid.gz

wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxcat.zip
unzip taxcat.zip

wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz

echo downloading and unzipping done

python $PWD/Diamond_BLAST_add_taxonomic_info/prepare_accession_to_description_db.py

echo four discription to accession number database done

#rm -rf nr.faa





