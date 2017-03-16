READme info for script to split fasta file into N files
==============================================

basic usage:

``./split_up_fasta_file.py`` -h 

$ split_up_fasta_file.py -n how_many_files -i in.fasta

this is a script to split up a big fasta file up into smaller ones.

this program will make a folder called split_fasta_files and
put the fasta files in there

why?

Can run many blast job instead of one
Some tools require 1 fasta entry per file


Options:
  -h, --help            show this help message and exit
  -n HOW_MANY_FILE, --num=HOW_MANY_FILE
                        How_many_file to split the fasta up into
  -i FILE, --in=FILE    in_file.fasta
