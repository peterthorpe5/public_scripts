READme: add taxonomic info and annotation to diamond blast ouptut
====================================================================

basic usage:

``./Diamond_blast_to_taxid.py`` -h 

warning: running to script uses a lot of RAM ~25GB. 

Main title: Diamond/BLAST-12 column -tab-blast to taxonomic-id info

can also work for 12 coloumn Blast output. - ALSO returns TOP-BLAST hits. Kingdom and Genus Distribution of these hits

purpose: script to pare the tabular output from $ diamond view -f tab -o name.tab
and get the description, tax id, species, kingdom information from NCBI taxonomy databse
why: diamond is SOOO fast. But does not include tax id info in the database.
author: Peter Thorpe September 2015. The James Hutton Insitute, Dundee, UK.
# NOTE: this will also work on standard blast output which does not have kingdom assignmnets.

Use as follows:

``./Diamond_blast_to_taxid.py -i diamond_tab_output -t /PATH_TO/NCBI_gi_taxid_prot.dmp -c /PATH/To/categories.dmp -n /PATH/To/names.dmp -d /PATH_TO_/description_database -o outfile.tab``

        or

``python Diamond_blast_to_taxid.py -i diamond_tab_output -p /PATH_TO/FILES -o outfile.tab``

Options:
  -h, --help            show this help message and exit
  -i FILE, --in=FILE    the tab output from diamond. use: $ diamond view -a
                        diamond.daa -f tab -o name.tab
  -p PATH, --path=PATH  Directory containing relevant taxonomy/database files
                        (set by -t, -c, -n, -d). Default is the current
                        workingdirectory. This is not used with the main input
                        and output filenames. This is where you put all the
                        downloadedNCBI files................... IF YOU GIVE
                        THE PATH YOU DONT NEED TO SET -t, -c, -n, -d)
  -t FILE, --taxid_prot=FILE
                        NCBI provided file access_tax_prot.dmp (from FTP site,
                        access_tax_prot.dmp.gz after unzipping). These file
                        required file options can be left blank if -p is
                        specified with a path to where all these can be found.
                        If -p /PATH/ is specified python will look in the
                        folder by default.
  -c FILE, --cat=FILE   NCBI provided kingdom catergories file categories.dmp
                        (from FTP site inside taxcat.zip).
  -n FILE, --names=FILE
                        NCBI provided names file names.dmp (from FTP site
                        inside taxdump.tar.gz
  -d FILE, --des=FILE   a databse of gi number-to-descrition. Generate either
                        using the shell script or by the following: export
                        BLASTDB=/PATH_TO/blast/ncbi/extracted # can only use
                        protein databases with DIAMOND. blastdbcmd -entry
                        'all' -db nr > nr.faa python
                        prepare_accession_to_description_db.py
  -o FILE, --out=FILE   Output filename - default=
                        infile_tab_blast_with_txid.tab


This script opens up a diamond tab output (-i) and looks up the relavant tax_id info (-t), look up the kingdom (-c)
looks for the species names (-n), looks up the descrtiption of the blast hit (-d):

FILES YOU NEED TO DOWNLOAD FOR THIS TO WORK!!
==============================================

	wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz.md5
	wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
	md5sum -c prot.accession2taxid.gz.md5
	gunzip prot.accession2taxid.gz

    Prot to tax_id: ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz).
    catergories file: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxcat.zip
    speices names: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

all files need to be uncompressed

do the following:
	wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz.md5
	wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
	md5sum -c prot.accession2taxid.gz.md5
	gunzip prot.accession2taxid.gz

    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxcat.zip
    unzip taxcat.zip

    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -zxvf taxdump.tar.gz

export BLASTDB=/path_to/blast/ncbi/

To generate the gi_to_des.tab databse:
blastdbcmd -entry 'all' -db nr > nr.faa

or download the nr database from, and the other 64 ish folders:
ftp://ftp.ncbi.nih.gov/blast/db/nr.00.tar.gz


YOU NEED TO PREPARE THE ACCESSION TO DESCRIPTION DATABASE:
==========================================================
``python prepare_accession_to_description_db.py --descriptions 2``

$ python prepare_accession_to_description_db.py -i nr.faa (default)-o acc_to_des.tab (dafault)

To prepare NCBI nr.faa:
blastdbcmd -entry 'all' -db nr > nr.faa

This script creates a database of accession number to description for the purpose
of post annotating diamond or BLAST tabular output.
Is it essential you prepare the databse before running the next script.

Use the --descriptions options to control how many descriptions for the same hit
are returned. default 4 hits.

descriptions = number of blast descriptions (default is 4) to return which are asscosiated with this accession number. 
More will be useful but your
file will get very big. Quickly!
./prepare_accession_to_description_db.py -h
Usage: Use as follows:

Options:
  -h, --help            show this help message and exit
  -i NR_FASTA_FILE, --fasta=NR_FASTA_FILE
                        nr_fasta_file, generate using blastdbcmd -entry 'all'
                        -db nr > nr.faa ,  you may need to: export
                        BLASTDB=/PATH/TO/ncbi/extracted default=nr.faa
  -n DESCRIPTIONS, --descriptions=DESCRIPTIONS
                        number of blast descriptions to return  which are
                        asscosiated with this GI number
  -o FILE, --out=FILE   Output filename: default=acc_to_des.tab


# Note: this return the top description in the NR discrition for each fasta file entry. 
# This can be modified to return "all". But the file will be much larger and therefore require more RAM


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
======================

By default this script will find the top hits by  1) Explicitly looking for the BLAST entry with the greatest bit score per query.
Script will also return the distribution of the kindgom and genus for these top hits.

Some notes on using Diamond:
=============================

# script to get the latest NR database and NT database and make a diamond blastdatabse.

# to install diamond from source
export BLASTDB=/PATH/TO/ncbi/extracted

``blastdbcmd -entry 'all' -db nr > nr.faa``

``/diamond-0.7.9/bin/diamond makedb --in nr.faa -d nr``

``diamond makedb --in uniprot_sprot.faa -d uniprot``

``diamond makedb --in uniref90.faa -d uniref90``

covert output to tab:
``diamond view -a diamond.daa -f tab -o name.tab``
