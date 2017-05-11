import os




count = 0



#f_in = open("R.padi_final_genome.v1.fasta", "r")

#f_out = open("R.padi_final_genome.v1_rename.fasta", "w")



f_in = open("hints2.rnaseq.ep.gff", "r")

f_out = open("hints3.rnaseq.ep.gff", "w")



for line in f_in:
    if "scaffol" in line:
        line = line.replace("|", "_")
        scaf,a,b,c,d,e,f,g,h = line.split("\t")
        scaf = scaf.split("_")[0]
        scaf = scaf.replace("scaffold", "RpS")
        data = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(scaf,a,b,c,d,e,f,g,h.rstrip())

        print >>f_out,data

