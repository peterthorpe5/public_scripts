#!/usr/bin/env python3
# title: reduce the info in the output of SNPeff
# the output is useful of course, but stops it being readble. 
# Author: Peter Thorpe. 


# imports
import sys
import os
from optparse import OptionParser
import datetime
import logging
import logging.handlers
import time


def gff(filename):
    "this is a function"
    old_to_new_dic = dict()
    old_tranID__prot = dict()
    expression_file = open(filename, "r")
    accession_to_new_name = dict()
    anno_out = open("V5_names_to_annotation.txt", "w")
    anno_out.write("#oldname\taccession\tgene_name\tannotation\n")
    seen_set = set([])
    for line in expression_file:
        id = line.split(";")[0]
        transID_line = line.split("transcript_id=gnl|WGS:PSQE|mrna.")[1]
        transID = transID_line.split(";product=")[0]
        prot_annot_temp = line.split(";product=")[1]
        prot_annot_temp = prot_annot_temp.split(";protein_id=")[0]
        prot_annot_full = prot_annot_temp
        prot_annot_temp = prot_annot_temp.replace("putative", "")
        prot_annot_temp = prot_annot_temp.replace("%", "")
        prot_annot_temp = prot_annot_temp.replace("/", "_")
        prot_annot_temp = prot_annot_temp.replace("(", "")
        prot_annot_temp = prot_annot_temp.replace(")", "")
        prot_annot_temp = prot_annot_temp.replace("-", "_")
        prot_annot = "_".join(prot_annot_temp.split())
        transID_prot = transID + "_" + prot_annot.rstrip()
        id = id.replace("ID=", "")
        old_gene = "gene" + id.replace("cds", "") 
        new_name = line.split("protein_id=")[1]
        new_name= new_name.rstrip()
        old_to_new_dic[old_gene.rstrip()] = new_name
        old_tranID__prot[old_gene.rstrip()] = transID_prot
        accession_to_new_name[new_name] = transID_prot
        if old_gene not in seen_set:
            anno_out.write("%s\t%s\t%s\t%s\n" % (old_gene, new_name, transID, prot_annot_full))
            seen_set.add(old_gene)
    expression_file.close()
    anno_out.close()
    #print(old_to_new_dic)
    #  'gene21027': 'RHN61598.1',
    return old_to_new_dic, old_tranID__prot, accession_to_new_name


def splitInfo(info, gene_prefix, old_to_new_dic):
    """func to reduc the info in the info coloumn
    I know this is reducing useful info, but people
    cannot digest that!!
    e.g. AC=2;AF=0.200;AN=10;DP=84;ExcessHet=1.5490;FS=0.000;MLEAC\
    =5;MLEAF=0.500;MQ=18.52;QD=25.57;SOR=4.615;\
    ANN=C|intergenic_region|MODIFIER|CHR_START-ACYPI\
    088670|CHR_START-ACYPI088670|intergenic_region|CHR_ST\
    ART-ACYPI088670|||n.39T>C||||||
    """
    modifier_type = info.split("|")[1]
    modifier = info.split("|")[2]
    geneextra = info.split("|")[3]
    gene = "NA"
    gene_prefix = "gene"
    try:
        gene = geneextra.split("-")[1]
    except:
        if geneextra == "":
            gene = "NA"
            return modifier_type, modifier, gene
        if geneextra.startswith(gene_prefix):
            gene = geneextra
        if len(geneextra.split(gene_prefix)) > 1:
            try:
                geneextra = geneextra.split(gene_prefix)[1]
                gene = gene_prefix + geneextra
            except:
                gene = "NA"
    return modifier_type, modifier, gene

def haplotype_reduce(incol):
    """reduce the info in the haplotype box"""
    try:
        haplotype = incol.split(":")[0]
    except:
        haplotype = ""
    return haplotype


def in_vcf(vcf, gene_prefix, old_to_new_dic, logger, outfile):
    """take in vcf file."""
    f_in = open(vcf, "r")
    f_out = open(outfile, "w")
    highImpact = outfile.split(".vcf")[0] + "_HIGH_IMPACT.vcf"
    f_highImpact = open(highImpact, "w")


    for line in f_in:
        if line.startswith("##contig"):
            continue
        if line.startswith("##scaffold"):
            continue
        if line.startswith("#"):
            if "INFO" in line:
                line = line.replace("INFO", "IMPACT")
                line = line.replace("FILTER", "SNP_location")
               
            f_out.write(line)
            f_highImpact.write(line)
            continue
        if not line.strip():
            continue # not a line anymore
        CHROM,  POS, ID,  REF, ALT, QUAL, FILTER,INFO, FORMAT, \
        RNIL_Control, SNIL_Control, RNIL_N116, SNIL_N116, RNIL_PS01, \
        SNIL_PS01= line.split("\t")
        modifier_type, modifier, gene = splitInfo(INFO, gene_prefix, old_to_new_dic)
        RNIL_Control = haplotype_reduce(RNIL_Control)
        SNIL_Control = haplotype_reduce(SNIL_Control)
        RNIL_N116 = haplotype_reduce(RNIL_N116)
        SNIL_N116 = haplotype_reduce(SNIL_N116)
        RNIL_PS01 = haplotype_reduce(RNIL_PS01)
        SNIL_PS01 = haplotype_reduce(SNIL_PS01)
        vcf_line = "\t".join([CHROM, POS, gene,
                             REF, ALT, QUAL,
                             modifier_type, modifier,
                             FORMAT, RNIL_Control,
                              SNIL_Control, RNIL_N116,
                              SNIL_N116, RNIL_PS01,
                              SNIL_PS01])
        f_out.write(vcf_line + "\n")
        if modifier == "HIGH":
            f_highImpact.write(vcf_line + "\n")
            
    f_out.close()
    f_in.close()
    f_highImpact.close()
        
        
        

usage = """ python reduce_vcf.py -h

python reduce_vcf.py --vcf snpef.vcf -o outreduced.vcf

reduces the info coloum to the genes only and 

# title: reduce the info in the output of SNPeff
# the output is useful of course, but stops it being readble. 
# Author: Peter Thorpe. 

"""



parser = OptionParser(usage=usage)


parser.add_option("--vcf", dest="vcf",
                  default="Medicago_V5_SNPs_freebayes.raw.g5mac3dp3.snp_effect.vcf",
                  help="vcf",
                  metavar="FILE")

parser.add_option("--logger", dest="logger",
                  default=None,
                  help="Output logger filename. Default: " +
                  "outfile_std.log",
                  metavar="FILE")

parser.add_option("--gene_prefix", dest="gene_prefix",
                  default="gen",
                  help="gene_prefix: eg. pea aphid is AC " +
                  ", helpful when looking for the genes")

parser.add_option("-o", "--out",
                  dest="outfile",
                  default="reduced.vcf",
                  help="Output filename (default: results.out)",
                  metavar="FILE")


(options, args) = parser.parse_args()



#--vcf
vcf = options.vcf
# gene_prefix
gene_prefix = options.gene_prefix
outfile = options.outfile


old_to_new_dic, old_tranID__prot, accession_to_new_name = gff("convert_name.txt")

#######################################################################
# Run as script
# Run as script
if __name__ == '__main__':
    # Set up logging
    if not options.logger:
        options.logger = "reduceSnpEff_vcf.log"
    logger = logging.getLogger('TranStart.py: %s' % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)
    try:
        logstream = open(options.logger, 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        outstr = "Could not open %s for logging" % options.logger
        logger.error(outstr)
        sys.exit(1)
    # Report input arguments
    logger.info(sys.version_info)
    logger.info("Command-line: %s", ' '.join(sys.argv))
    logger.info("Starting testing: %s", time.asctime())
    if not os.path.isfile(vcf):
        print("vcf not here motherfucker!")
    in_vcf(vcf, gene_prefix, old_to_new_dic, logger, options.outfile)
