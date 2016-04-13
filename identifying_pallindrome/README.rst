READme for finding pallindromes
==============================================

Requirements:
Biopython, numpy, scipy, reportlab(if you want to draw the genome diagrams)


Basic usage for all scripts. 
===========================
``./name_of_script.py`` -h 


This has multiple step:

step 1: Indefity pallindromes
============================

``./identify_pallidromes.py`` --fasta /path/file.fasta --folder (where_ .fasta is default= os.getcwd()) --base_name (species) --min_pallindrome_len (default 7) --max_pallindrome_len (default 1000)

e.g.  ``./identify_pallidromes.py`` --fasta Gr4 --folder /home/nematode/Gr4 --base_name Gr4

This is a master script to call the palindrome library. This will identify
pallindromes the come under a certain probability of occuring at random and
allows missmtaches...

The search will look for pallindromes from min_length to max_length

The function called are found in
``./palindrome_library.py`` 


step 2: reduce pallindrome redundancy
====================================

``./pallindrome_repeat_reducer.py`` --folder /path/to/ (default= os.getcwd()) --base_name (name_for_species)

e.g. ``./pallindrome_repeat_reducer.py`` --base_name Gr4 --folder /home/misc_python/identifying_pallindrome_final/Gr4

  TITLE: this is a function to reduce the number of repeated palindromes found in a txt file. Repeated palindromes are: Palindromes that 'fit'
  inside another palindrome sequence. For example p and m. m is the longer sequence, if p start > m start and p end <  m end, p fits in m. From our
 results file we need to reduce palindromes that 'fit' inside another...

 
 
 step 3: pallindromes: genic_or_non-genic?
 ==============================================
 
 ``./genic_or_non-genic.py`` --ptt  --gff --folder --base_name
 e.g. ``./genic_or_non-genic.py`` --folder /home/nematode/mitochondria/Gr1 --base_name Gr1 --gff /home/nematode/mitochondria/genbank_files/Gr1.gff
 
  e.g. ``./genic_or_non-genic.py`` --gff Gr4.gff --base_name Gr4

 #   TITLE: This is a funtion to compare the results of the indentified palindromes
#against the starts and end positions of the genome in question.

#why?: if the palindrome occurs within the boundries of a gene start....end
#position, then this is said to be a genic palindrome. Palindromes are expected
#to be found in intergenic regions, not within the boundires of the gene.
#The result from this code will help us with the hypothesis, with the comparison



step 4: pallindrom_clustering_assessment.py
==========================================
 ``./pallindrom_clustering_assessment.py`` -h
 
 e.g.  ``./pallindrom_clustering_assessment.py`` --folder /home/nematode/mitochondria/Gr1 --base_name Gr1
 
 e.g. ``../pallindrom_clustering_assessment.py`` --base_name Gr4

 this is a function to return start position minus the next start position
using an excisting data set.

WHY?: This will provide us with the data to see if there is any evidence of
clustering of palindromes against a random data set.I would predict the 'random'
to have a normal distribution of distance between one palindrome to the
next.... if we compare our real data against the data generated from the
random data we can prove if palindrome clustering is happening.

#   TITLE: This is a function to return start position minus start positions
#in an ordered data set of identified palindromes. This is fo the identification
#of palindrome clustering



step 5:
=======
do some stats on the output!


step 6: drawing_with_palindromes.py
===================================
# requires report lab

TO DO modify script. 