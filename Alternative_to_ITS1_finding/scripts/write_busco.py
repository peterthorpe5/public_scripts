import os

for filename in os.listdir(".") :
    if not filename.endswith(".fna") : continue
    name_prefix = filename.split("_genomic")[0]
    name_out = "busco_%s.sh" %(name_prefix)
    f_out= open(name_out, 'w')

    cd = '#$ -l hostname="n13*"\n#mkdir LINEAGE\n#ln -s ~/scratch/Downloads/BUSCO_v1.1b1/Eukaryota ./LINEAGE/\ncd /home/pt40963/scratch/tree_health/Phytophthora_genomes\n'
    f_out.write("%s" % (cd))

    command = "python3 /home/pt40963/scratch/Downloads/BUSCO_v1.1b1/BUSCO_v1.1b1.py -in %s -l ./LINEAGE/Eukaryota -o busco_%s  -m genome -f -Z 827000000 --long --cpu 4 " % (filename, name_prefix)

    f_out.write(command +"\n")
    f_out.close()
print "Done"
