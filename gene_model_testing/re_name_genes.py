#!/usr/bin/env python
# Title:
# script to reformt gene names after gene prediction.
# GT sometime returns . as a name for all genes.
import sys
from optparse import OptionParser


def open_file(in_file):
    """funtion to open the tab sep file.
    returns a \n separeted list
    read in the tab file. Reads whole file into memory.
    Could be altered for more efficiency
    """
    with open(in_file) as file:
        return file.read().split("\n")




def reformat_genome_tools_gtf_column(gtf, out_prefix):
    """function to rename the gene in a GTF file. The automatically
    gnerated Augustus names/Braker name need to be altered.
    These are altered to perfix00001 - then 5 interger number gene name"""
    f_in = open_file(gtf)
    outfile = out_prefix + ".gtf"
    f_out = open(outfile, 'w')
    gene_count = 0
    for line in f_in:
        if line.startswith("#"):
            f_out.write(line.rstrip() + "\n")
            continue
        # split the coloumn up in the gff based on \t
        scaf, aug, info, start, stop, stats, direction,\
                    more_info, gene_info = line.split("\t")
        #count the genes
        if info == "gene":
            gene_count = gene_count + 1
        # new gene names to 5 interger places.
        # a long winded way of getting to the original g number.
        # this has to be long winded as the line in the GFF do not follow
        # one rule.
        original_gene_num = gene_info.split('"')[1]
        original_gene_num = original_gene_num.split(";")[0]
        original_gene_num = original_gene_num.split('g')[1]
        #temp variable to the: =g1 or whatever the current gene is
        current_gff_gene = "g" + original_gene_num
        original_gene_num = int(original_gene_num)

        temp_gene_name = gene_info.replace(current_gff_gene, "=g")

        gene_number = "=%s%05d" %(prefix, original_gene_num)
        #if info == "gene":
            #print gene_number
        new_gene_names = temp_gene_name.replace("=g", gene_number)
        new_gene_names = new_gene_names.replace("=", "")

        data = "\t".join([scaf,
                          aug,
                          info,
                          start,
                          stop,stats,\
                          direction,
                          more_info,\
                          new_gene_names + "\n"])
        f_out.write(data)
    f_out.close()
    f_in.close()


def reformat_gtf_column(gtf,out_prefix):
    """finction to rename the gene in a GTF file. The automatically
    gnerated Augustus names/Braker name need to be altered.
    These are altered to perfix00001 - then 5 interger number gene name"""
    f_in = open(gtf, "r")
    outfile = out_prefix+".gff"
    f_out = open(outfile, 'w')
    gene_count = 0
    for line in f_in:
        if line.startswith("#"):
            print >>f_out,line.rstrip()
            continue
        # split the coloumn up in the gff based on \t
        scaf,aug,info,start,stop,stats,direction,\
                            more_infor,gene_info = line.split("\t")
        #count the genes
        if info == "gene":
            gene_count = gene_count+1

        # new gene names to 5 interger places.
        # a long winded way of getting to the original g number.
        # this has to be long winded as the line in the GFF do not follow
        # one rule.
        original_gene_num = gene_info.split(".t")[0]
        original_gene_num = original_gene_num.split(";")[0]
        original_gene_num = original_gene_num.split('g')[1]

        #temp variable to the: =g1 or whatever the current gene is
        current_gff_gene = "g"+original_gene_num

        #sanity check
        assert int(original_gene_num.rstrip()) == gene_count

        temp_gene_name = gene_info.replace(current_gff_gene, "=g")
        gene_number = "=%s%05d" %(prefix, gene_count)
        #if info == "gene":
            #print gene_number
        new_gene_names = temp_gene_name.replace("=g", gene_number)
        new_gene_names = new_gene_names.replace("=", "")


        data = "\t".join([scaf,
                          aug,
                          info,
                          start,
                          stop,stats,\
                          direction,
                          more_info,\
                          new_gene_names + "\n"])
        f_out.write(data)
    f_out.close()
    f_in.close()


def reformat_gff_column(gff, prefix, out_prefix):
    """function to rename the gene in a GFF file. The automatically
    gnerated Augustus names/Braker name need to be altered.
    These are altered to perfix00001 - then 5 interger number gene name"""
    f_in = open(gff, "r")
    outfile = out_prefix + ".gff"
    f_out = open(outfile, 'w')
    gene_count = 0
    for line in f_in:
        if line.startswith("#"):
            f_out.write(line.rstrip() + "\n")
            continue
        # split the coloumn up in the gff based on \t
        scaf, aug, info, start, stop, stats, direction,\
                    more_info, gene_info = line.split("\t")
        #count the genes
        if info == "gene":
            gene_count = gene_count + 1
            gene_number = "ID=%s%05d" %(prefix, gene_count)
        gene_number = "%s%05d" %(prefix, gene_count)
        data = "\t".join([scaf,
                          aug,
                          info,
                          start,
                          stop,stats,\
                          direction,
                          more_info,\
                          gene_number + "\n"])
        f_out.write(data)
    f_out.close()
    f_in.close()



def reformat_fasta_gene_name(filename, prefix, out_prefix):
    """this function re-write a file as a fasta file but with
    altered names.  the default g. in augustus output is replace by
    prefix. Gene numbers are to
    5 intiger places. """
    outfile = out_prefix + "_alt.fasta"
    f= open(outfile, 'w')
    f_in = open(fasta, "r")
    prefix = str(prefix)
    count = 0
    for seq_record in SeqIO.parse(filename, "fasta"):
        print(count)
        count = count + 1
        gene_name_id = str(seq_record.id).split("|")[0]
        gene_name_id = gene_name_id.split("gene=")[0]
        # magic to rename the gene
        gene_number = "%s%05d" %(prefix, count)
        gene_name = gene_name_id.replace("g", gene_number)
        seq_record.id = gene_name
        seq_record.name = ""
        seq_record.description = ""
        SeqIO.write(seq_record, f, "fasta")
    f.close()
    f_in.close()


def reformat_hints_scaffold_name(hints, prefix, out_prefix):
    "function to reformat name. Remove pipes and make names shorter"
    f_in = open(hints, "r")
    outfile = out_prefix+"_alt.hints"
    f_out = open(outfile, 'w')
    prefix = str(prefix)
    for line in f_in:
        if "scaffol" in line:
            #line = line.replace("_", "_")
            scaf,a,b,c,d,e,f,g,h = line.split("\t")
            scaf = scaf.split("|")[0]
            scaf = scaf.replace("scaffold", prefix)
            data = "\t".join([scaf,
                              aug,
                              info,
                              start,
                              stop,stats,\
                              direction,
                              more_info,\
                              gene_number + "\n"])
            f_out.write(data)
    f_out.close()
    f_in.close()


#################################################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

$ python re_format_fasta_hints -f genome.fasta --hints hints_file --prefix Mc

script to reformt scaffold names for Braker. It doesnt seems to like '|' or long names

requires Biopython
"""

parser = OptionParser(usage=usage)

parser.add_option("-f","--fasta", dest="fasta", default=False,
                  help="fasta file to have names altered")
parser.add_option("--hints", dest="hints", default=False,
                  help="hintsfile",
                  metavar="FILE")
parser.add_option("--gff", dest="gff", default=False,
                  help="gff file",
                  metavar="FILE")
parser.add_option("--gtf", dest="gtf", default=False,
                  help="gtf file",
                  metavar="FILE")
parser.add_option("--gtf_gt", dest="gtf_gt", default=False,
                  help="gtf made by genome tools file",
                  metavar="FILE")
parser.add_option("-p", "--prefix", dest="prefix", default=None,
                  help="prefix to alter the gene names to "
                  " e.g. Rpa00001",
                  metavar="FILE")
parser.add_option("-o", "--out_prefix", dest="out_prefix", default=None,
                  help="prefix to the output filenames")


(options, args) = parser.parse_args()

fasta = options.fasta
hints = options.hints
gff = options.gff
prefix = options.prefix
out_prefix = options.out_prefix
gtf = options.gtf
gtf_gt = options.gtf_gt




# Run as script
if __name__ == '__main__':
    if gtf_gt:
        reformat_genome_tools_gtf_column(gtf_gt, out_prefix)

    #run the program

    #biopython imports
    if hints:
        reformat_hints_scaffold_name(hints, prefix, out_prefix)
    if gff:
        reformat_gff_column(gff, prefix, out_prefix)

    if gtf:
        reformat_gtf_column(gtf,out_prefix)

    if fasta:
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        reformat_fasta_gene_name(fasta, prefix, out_prefix)

