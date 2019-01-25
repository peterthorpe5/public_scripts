import os
from sys import stdin,argv


def seq_getter(filename_in, wantedfile, outfile):
    "this is a function to open up a .xml file blast results, the out put of\
is all the unique hit"
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    f= open(outfile, 'w')
    wanted = open(wantedfile, "r")

    names = wanted.readlines()
    #print names
    wanted_data = [line.replace("\t", "").rstrip("\n") for line in names
              if line.strip() != ""]
    name_set = set([])
    for i in wanted_data:
        if not i.startswith("#"):
            name_set.add(i)
    #print wanted_data

    cds_database = SeqIO.index(filename_in, "fasta")
    #record = SeqIO.read(filename, "fasta")
    for i in name_set:
        if "\r\n" in i:
            i = i.replace("\r\n","")
        #print i
        seq_record = cds_database[i]
        #print 'boomshanka'
        SeqIO.write(seq_record, f, "fasta")
    f.close()
    return True


seq_getter(argv[1],argv[2], argv[3])

#seq_getter('assembly2_scaffolds.fasta',\
           #'scaffold318.fasta')
print 'done'

