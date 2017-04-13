READme THIS IS FOR LEGACY NCBI FORMAT - THE LATEST DATA WILL NOT WORK WITH THESE SCRIPTS. BUT OLDER DATA WILL ONLY WORK WITH THESE SCRIPT
=========================================================================================================================================

basic usage:

``./Diamond_blast_to_taxid.py`` -h 

warning: running to script uses a lot of RAM ~25GB. 

Main title: Diamond/BLAST-12 column -tab-blast to taxonomic-id info

can also work for 12 coloumn Blast output. - ALSO returns TOP-BLAST hits. Kingdom and Genus Distribution of these hits

purpose: script to pare the tabular output from $ diamond view -f tab -o name.tab
and get the description, tax id, species, kingdom information from NCBI taxonomy databse

why: diamond is SOOO fast. But does not include tax id info in the database.

author: Peter Thorpe September 2015. The James Hutton Insitute, Dundee, UK.


Use as follows:

``./Diamond_blast_to_taxid.py -i diamond_tab_output -t /PATH_TO/NCBI_gi_taxid_prot.dmp -c /PATH/To/categories.dmp -n /PATH/To/names.dmp -d /PATH_TO_/description_database -o outfile.tab``

        or

``python Diamond_blast_to_taxid.py -i diamond_tab_output -p /PATH_TO/FILES -o outfile.tab``


This script opens up a diamond tab output (-i) and looks up the relavant tax_id info (-t), look up the kingdom (-c)
looks for the species names (-n), looks up the descrtiption of the blast hit (-d):

# NOTE: this will also work on standard blast output which does not have kingdom assignmnets.

    Prot to tax_id: ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz).
    catergories file: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxcat.zip
    speices names: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

all files need to be uncompressed

do the following:

    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz
    gunzip gi_taxid_prot.dmp.gz

    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxcat.zip
    unzip taxcat.zip

    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -zxvf taxdump.tar.gz


To generate the gi_to_des.tab databse:
blastdbcmd -entry 'all' -db nr > nr.faa


python ~/misc_python/diamond_blast_to_kingdom/prepare_gi_to_description_databse.py

# Note: this return the top description in the NR discrition for each fasta file entry. This can be modified to return "all". But the file will be much larger and therefore require more RAM


    BLAST DATA we be returned as:

qseqid = Query Seq-id (ID of your sequence)
sseqid = Subject Seq-id (ID of the database hit)
pident = Percentage of identical matches
length = Alignment length
mismatch = Number of mismatches
gapopen = Number of gap openings
qstart = Start of alignment in query
qend = End of alignment in query
sstart = Start of alignment in subject (database hit)
send = End of alignment in subject (database hit)
evalue = Expectation value (E-value)
bitscore = Bit score
salltitles = TOP description of the blast hit
staxids = tax_id
scientific_name
scomnames = common_name
sskingdoms = kingdom


TOP BLAST HITS FINDER:
By default this script will find the top hits by two methods. 1) assuming order in BLAST out file 2) Explicitly looking for the BLAST entry with the greatest bit score per query.
Script will also return the distribution of the kindgom and genus for these top hits.



Some notes on using Diamond:


# script to get the latest NR database and NT database and make a diamond blastdatabse.


# to install diamond from source
export BLASTDB=/PATH/TO/ncbi/extracted


``blastdbcmd -entry 'all' -db nr > nr.faa``

``/diamond-0.7.9/bin/diamond makedb --in nr.faa -d nr``

``diamond makedb --in uniprot_sprot.faa -d uniprot``

``diamond makedb --in uniref90.faa -d uniref90``

covert output to tab:
``diamond view -a diamond.daa -f tab -o name.tab``


