import os
from sys import stdin,argv


def seq_getter(filename_in, wantedfile, threshold, outfile):
    "this is a function to open up a .xml file blast results, the out put of\
is all the unique hit"
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    f= open(outfile, 'w')
    wanted = open(wantedfile, "r")
    #len of seq to filter at
    threshold = int(threshold)

    names = wanted.readlines()
    #print names
    wanted_data = [line.replace("\t", "").rstrip("\n") for line in names
              if line.strip() != ""]
    name_set = set([])
    for i in wanted_data:
        if not i.startswith("#"):
            name_set.add(i)
    #print wanted_data
    #record = SeqIO.read(filename, "fasta")
    for seq_record in SeqIO.parse(filename_in, "fasta"):
        #print record.description
        
        if not seq_record.id in name_set:
            if len(seq_record.seq)> threshold:
                SeqIO.write(seq_record, f, "fasta")
    f.close()
    return True


seq_getter(argv[1],argv[2], argv[3], argv[4])

#seq_getter('assembly2_scaffolds.fasta',\
           #'scaffold318.fasta')
print 'done'

