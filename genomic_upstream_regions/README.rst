READme info for script to get upstream regions
==============================================

basic usage:

``./get_upstream_regions.py`` -h 

script to get the upstream regions of genes of interest 
#script will return upt to the gene if the full length falls within that gene.
#also, script will retunr reverse complemnet of negative strand coded genes.

Usage: Use as follows:

$ python get_upstream_regions.py --coordinates
        coordinate_file.fasta -g genome_sequence.fasta
        -upstream <int> number of nucleotides upstream of strat of gene to return
        e.g.  -u 1000
        -z user_defined_genic (how much of the gene to return)
        -o outfile_name

Requirements:
python 2.7 and biopyton.  
 (print >> is 2.X only).

This will return (--upstream number) of nucleotides to the start of your genes(s) of
interest (-g) gene_file using data from (-c). Gene file can either be space, tab or
 separated..

The coordinate file can be generated using a GFF3 file and a linux command line using:

grep "gene" name.gff3 | grep -v "#" | cut -f1,4,5,7,9 > format_for_py_script.out


yeilding this resulting file:

scaffold        start   stop    strand(+/-)     ID=gene_number
GROS_00001      2195    3076    -       ID=GROS_g00002
GROS_00001      8583    10515   +       ID=GROS_g00005.....

The script will check that the start is always less than the end. GFF file should have
starts < stop irrespective of the coding direction

To get all the genes in the file do:

cut -f5 format_for_py_script.out > all_gene_names.out

MORE help / example:

This is an example I ran for G. pallida:

python get_upstream_regions.py -c format_for_py_script.out -g Gpal.v1.0.fas
-f all_gene_names.out -u 225 -z -125
-o Gp.all_gene_names.out_225up_125genic.fasta
> warning_all_gene_names.out_225up_125genic.out


This got (-u 225) 225 bp upstream and (-z 125) 125bp into the current gene for all the
genes in (-f) all_gene_names.out. By default -z is zero. So you dont need to specify this,
unless you specifically want a piece of the current gene being searched for.


Options:
  -h, --help            show this help message and exit
  -c COORDINATE_FILE, --coordinates=COORDINATE_FILE
                        NOTE: coordinate_file can generate using linux command
                        line of GFF file:  grep 'gene' name.gff3 | grep -v '#'
                        |  cut -f1,4,5,7,9 > format_for_py_script.out .Default
                        = format_for_py_script.out
  -g GENOME_SEQUENCE, --genome=GENOME_SEQUENCE
                        genome_sequence.fasta  -  this has to be the file used
                        to generate the gene models/GFF file
  -f GENES_FILE, --gene_names=GENES_FILE
                        a file with a list of gene names to get the upstream
                        regions for
  -u UPSTREAM, --upstream=UPSTREAM
                        the amount of nucleotide upstream of the gene start,
                        taking into account gene directions, to return in the
                        outfileby default this will not return sequences of
                        50bp or less. If you require these alterlines 204 and
                        214
  -z USER_DEFINED_GENIC, --user_defined_genic=USER_DEFINED_GENIC
                        the number of nucleotides from within the gene to
                        return, default is 0
  -o FILE, --output=FILE
                        Output filename (fasta file)


						
``./get_upstream_regions.py`` --coordinates coordinate_file.fasta -g genome_sequence.fasta 
-upstream <int> number of nucleotides upstream of strat of gene to return 
e.g.  -u 1000
-z user_defined_genic (how much of the gene to return) -o outfile_name

This will return (--upstream number) of nucleotides to the start of your genes(s) of interest 
(-g) gene_file using data from (-c). Gene file can either be space, tab or \n separated..

The coordinate file can be generated using a GFF3 file and a linux command line using:

``grep`` "gene" name.gff3 | cut -f 1,4,5,7,9 > format_for_py_script.out
WARNING: the grep method can return lots of unwanted lines. Use with caution

or use the python script:
``./re_format_GFF_Mcscanx.py`` -h re_format_GFF_Mcscanx.py


yeilding this reulting file:
scaffold	start	stop	strand(+/-)	ID=gene_number
GROS_00001	2195	3076	-	ID=GROS_g00002
GROS_00001	8583	10515	+	ID=GROS_g00005.....

The script will check that the start is always less than the end. GFF file should have starts < stop irrespective of the coding direction



	MORE help / example:

This is an example I ran for G. pallida:

``python get_upstream_regions.py`` -c format_for_py_script.out -g Gpal.v1.0.fas -f all_gene_names.out 
-u 225 -z -125 -o Gp.all_gene_names.out_225up_125genic.fasta 
> warning_all_gene_names.out_225up_125genic.out


This got (-u 225) 225 bp upstream and (-z 125) 125bp into the current gene for all the genes in (
-f) all_gene_names.out. By default -z is zero. So you dont need to specify this, unless you specifically 
want a piece of the current gene being searched for.

