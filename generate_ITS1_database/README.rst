READme Generate ITS1 database
=============================

basic usage:

STEPS
=====

1) Open the ``excel sheet`` and save the file as a ``txt file``, e.g. ``Phytophthora_ITS1_DB_v0.006_20180620.txt``.


2) You have to update the following in the shell to run it:  open ``./Auto_make_database.sh`` in a text editor, put the name 
in DATA_TAB_FILE variable, WITHOUT a space, and the version you want it to be.  Eg. 

# PUT the db file name here

``DATA_TAB_FILE="Phytophthora_ITS1_DB_v0.006_20180620.txt"``

``version="0.006"``

ALTER this line:

``source ~/misc_python/THAPBI/THAPBI-pycits/venv-THAPBI-pycits-py3/bin/activate``

This is my VM. Put a pth to yours, if you want to use a VM


3) Run the shell on Linux:  ``./Auto_make_database.sh``

each python script in ./bin/ can be run with -h to see what they do



Requires:
=========
python3

Biopython

hmmsearch

muscle

internet access

OUTPUT
======
The main fasta file on the ITS1 only region called: 

``Phytophora_ITS_database_v${version}.fasta``


``seqid_domain_coordinates.txt`` = the coordinates of the ITS1 domain as Identified by the HMM.

Folder -  database_files
========================
``Duplicate.WARNINGS`` = This contains VERY imporatant warning files. Such as ducplicated seqeucnes, which have been combined into
one sequence. Sequences which have illegal characters, these are rejected by the program and NOT put in the 
final fasta file. 

``Duplicate2.WARNINGS`` = This is after the second round of filtering of duplicated sequences and checking for
Illegal charaters. This should be zero for all enteries.

``check_me_for_completeITS_may_be_a_lie.txt`` = sequences identified as ITS1, but not called by HMM after filtering. 
This sometimes happens if the domain is the whole sequence. You may have to manually update the DB, if you think this is real. 

``ITS1_hmm_domain_table.out`` -  HMM output from the ITS1 search. Why? The program downloads whatever seq is availble for that accession.
This may be a full region ~700nt. So We have to identify the ITS1 region. 

``Phytophthora_ITS1_DB_v${version}_full_len.WARNING`` - This check the name you have in your txtfile/Excel database an compares it to the name
 from GenBank. Please check this make sure accssion numbers have not been given to the wrong species by mistake. This happened once.
 This file should help identifiy these problems. 
 
``Phytophthora_ITS1_DB_v${version}_full_len.fasta`` - This is the full seq in the Gnebank record for the accession of interest. 

``Phytophthora_ITS1_DB_v0.006_20180620_ITS_only.fasta.aln`` - This is an alignmnet of the ITS1 region where the user can draw a tree of it.

Version 0.003 change.
=====================
20 base pair off each sequence removed. This si the right primer sequence which is in the 5.8 region.
P. austr spelling mistake changed. 

Version 0.004
=============
P. alni update

version 0.005
=============
P. gonap updated

Version 0.006
=============
updated by David cooke. 
Phytophthora castanetorum
Phytophthora tubulina
Phytophthora tyrrhenica
Phytophthora vulcanica
Phytophthora gregata
Phytophthora oleae
Phytophthora condilina
Phytophthora balyanboodja
Phytophthora pseudorosacearum
Phytophthora kwongonina 
Phytophthora cooljarloo
Phytophthora versiformis
