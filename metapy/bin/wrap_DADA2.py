#!/usr/bin/env python3
#
# title: wrap and run DADA2 pipeline
# https://benjjneb.github.io/dada2/tutorial.html
# requires R and DAD2
# in the cluster
# author: Peter Thorpe, Leighton Pritchard September 2017.
# The James Hutton Insitute, Dundee, UK.
# imports
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import argparse
import subprocess
import sys
import gzip
from Bio import SeqIO


####################################################################
# to run the script

usage = """usage :

run DAD2 pipline with metapy.py
python wrap_DAD2py -l left.fq.gz -r right.fq.gz

requires:
 R and DAD2

 source("https://bioconductor.org/biocLite.R")
biocLite("dada2")

"""

if sys.version_info[:1] != (3,):
    # e.g. sys.version_info(major=3, minor=6, micro=7,
    # releaselevel='final', serial=0)
    # break the program
    print ("currently using:", sys.version_info,
           "  version of python")
    raise ImportError("Python 3.x is now required for LTG.py")
    print ("did you activate the virtual environment?")
    print ("this is to deal with module imports")
    sys.exit(1)

VERSION = "Pycits classify OTU: v0.0.2"
if "--version" in sys.argv:
    print(VERSION)
    sys.exit(1)


def make_folder(WORKING_DIR, folder, exist_ok=True):
    """function to make a folder with desired name"""
    dest_dir = os.path.join(WORKING_DIR, folder)
    try:
        os.makedirs(dest_dir)
    except OSError:
        print ("folder already exists " +
               "I will write over what is in there!!")
    return dest_dir

def convert_file_to_tuple(infile):
    """factor out old function to parse and convert
    R seq table"""
    f_in = open(infile, "r")
    seqs = []
    abundance = []
    count = 1
    lines = f_in.readlines()
    seq_line = lines[0]
    seq_line = seq_line.replace('"', '').rstrip()
    seq_line = seq_line.split(" ")
    for sequence in seq_line:
        seqs.append(sequence)
    abundances = lines[1]
    abundances = abundances.split(" ")
    # first element if not wanted, as this is the read name.
    for abunda in abundances[1:]:
        abundance.append(abunda.rstrip())
    abund_seq = zip(abundance, seqs)
    return abund_seq


def write_as_fasta(infile, right_trim, outfile):
    """func to open up the seq table and write out as fasta
    File is formated as two lines.
    line1 = a seris of seq string.
    line2 = DNAMIX 942 488 336 28 """
    abund_seq = convert_file_to_tuple(infile)
    f_out = open(outfile, "w")
    count = 0
    for entry in abund_seq:
        count = count + 1
        abundance = entry[0]
        seq = entry[1]
        name = "seq%d_%s" % (count, abundance)
        seq_record = SeqRecord(Seq(seq[:-right_trim]),
                               id=name, name="",
                   description="")
        SeqIO.write(seq_record, f_out, "fasta")


def decompress(infile):
    """function to decompress gzipped reads"""
    cmd = ' '.join(["gunzip", infile])
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    return os.path.splitext(infile)[0]


def compress(infile):
    """function to compress reads, make them .gz"""
    cmd = ' '.join(["gzip", "-f", infile])
    pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)


def filter_my_fastq_file (in_fastq, trim_len, out_fastq):
    """function to parse fastq files and trim off a set
    lenght (trim_len)
    Writes to a new fq file"""
    # creat a new fastq file to write to
    out_file = open(out_fastq, "w")
    # open the fastq file
    in_file = open(in_fastq)
    # iterate through the fastq file
    for title, seq, qual in FastqGeneralIterator(in_file):
        out_file.write("@%s\n%s\n+\n%s\n" % (title,
                                             seq[trim_len:],
                                             qual[trim_len:]))
    out_file.close()
    in_file.close()



def run_sub_pro(command):
    """function to run subprocess as this will be repeated
    Takes in a command
    """
    pipe = subprocess.run(command, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)


def get_args():
    parser = argparse.ArgumentParser(description="wrap DAD2 " +
                                     " clusters :  %s " % usage,
                                     add_help=False)
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-l", "--left", dest="left",
                          default=None,
                          help="left read file")


    optional.add_argument("-r", "--right", dest="right",
                          default=None,
                          help="right read file")

    optional.add_argument("-w", "--wd", dest="working_dir",
                          default=os.getcwd(),
                          help="working_dir")

    optional.add_argument("-d", "--db", dest="database",
                          default=None,
                          help="database")

    optional.add_argument("--left_trim", dest='left_trim',
                          action="store",
                          default=53,
                          type=int,
                          help="left_trim for primers or conserved " +
                          "regions. Default 53 ")

    optional.add_argument("--right_trim", dest='right_trim',
                          action="store",
                          default=20,
                          type=int,
                          help="right_trim for primers or conserved " +
                          "regions.  Default 20")

    optional.add_argument("--compress", dest='compress',
                          action="store",
                          default="No",
                          type=str,
                          help="compress the read files? " +
                          " Default No")

    optional.add_argument("-h", "--help",
                          action="help",
                          default=argparse.SUPPRESS,
                          help="Displays this help message"
                          " type --version for version")
    args = parser.parse_args()
    return args

args = get_args()
left_trim = args.left_trim
right_trim = args.right_trim
old_right = os.path.abspath(args.right)
old_left = os.path.abspath(args.left)
working_dir = args.working_dir
database = args.database
READ_PREFIX = os.path.split(old_right)[-1].split("_R")[0]
COMPRESS = args.compress

# run the program
if __name__ == '__main__':
    # remove the primers from the fq files
    read_directory = os.path.split(old_left)[:-1][0]
    # decompress these for biopython
    if old_left.endswith(".gz"):
        old_left = decompress(old_left)
        old_left = old_left.split(".gz")[0]
    if old_right.endswith(".gz"):
        old_right = decompress(old_right)
        old_right = old_right.split(".gz")[0]
    # make the filtered fq files folder
    file_directory = os.path.split(old_left)[:-1][0]
    make_folder(file_directory, "filtered")
    # remove primers
    left = "%s_primers_trimmed.fastq" % (os.path.split(old_left)[-1].split(".fastq")[0])
    right = "%s_primers_trimmed.fastq" % (os.path.split(old_right)[-1].split(".fastq")[0])
    filter_my_fastq_file(old_left, left_trim,
                         os.path.join(read_directory,
                                      left))
    filter_my_fastq_file(old_right, right_trim,
                         os.path.join(read_directory,
                                      right))
    # compress these files.
    # TO DO use python gzip lib
    # uncompressed = [os.path.join(read_directory, left),
                    # os.path.join(read_directory, right)]
    compress(os.path.join(read_directory,
                          left))
    compress(os.path.join(read_directory,
                          right))
    # Set up Rscript file
    shell = open("dada2.R", "w")
    make_folder(working_dir, "dada2")
    shell.write("#!Rscript\n")
    shell.write('library(dada2); packageVersion("dada2")\n')
    print("read_directory = ", read_directory)
    r_path = 'path <- "%s"' % (read_directory)
    shell.write(r_path)
    shell.write("\n")
    shell.write('filtpath <- file.path(path, "filtered")\n')
    shell.write('fns <- list.files(path, pattern="primers_trimmed.fastq.gz")\n')
    # TODO "maxEE=3, truncQ=2 -  put these as command line options
    f_and_t = " ".join(["filterAndTrim(file.path(path,fns),",
                        "file.path(filtpath,fns),",
                        "maxEE=2, truncQ=2, rm.phix=TRUE,",
                        "compress=TRUE, verbose=TRUE, multithread=TRUE)\n"])
    shell.write(f_and_t)
    # CHANGE if different file extensions
    filt_path = os.path.join(file_directory, "filtered")
    fil_out = 'filtpath <- "%s"' % filt_path
    shell.write(fil_out)
    shell.write("\n")
    # Assumes filename = sample_XXX.fastq.gz
    shell.write('filts <- list.files(filtpath, pattern="primers_trimmed.fastq.gz", full.names=TRUE)\n')
    shell.write('sample.names <- sapply(strsplit(basename(filts), "_"), `[`, 1)\n')
    shell.write('names(filts) <- sample.names')
    shell.write("\n")
    # Learn error rates
    shell.write('set.seed(100)')
    shell.write("\n")
    shell.write("err <- learnErrors(filts, multithread=TRUE, randomize=TRUE)\n")
    shell.write("plotErrors(err, nominalQ=TRUE)\n")
    # Infer sequence variants
    shell.write('dds <- vector("list", length(sample.names))\n')
    shell.write('names(dds) <- sample.names\n')
    shell.write('for(sam in sample.names) {\n')
    ment = "n"
    temp = '\tcat("Processing:", sam, "\%s")' % (ment)
    shell.write(temp + "\n")
    shell.write('\tderep <- derepFastq(filts[[sam]])\n')
    shell.write('\tdds[[sam]] <- dada(derep, err=err, multithread=TRUE)\n}\n')
    # Construct sequence table and write to disk
    shell.write('seqtab <- makeSequenceTable(dds[sam])\n')
    seqtabs_dir = os.path.join(working_dir, "dada2", "seqtab.rds")
    seqtabs_dir_out = 'saveRDS(seqtab, "%s")\n' % seqtabs_dir
    shell.write(seqtabs_dir_out)
    read_in = 'st1 <- readRDS("%s")\n' % seqtabs_dir
    shell.write(read_in)
    # Remove chimeras
    chmira = 'st1, method="pooled", multithread=TRUE'
    chmira_out = 'seqtab <- removeBimeraDenovo(%s)' % chmira
    shell.write(chmira_out)
    shell.write("\n")
    seq_out = 'write.table(seqtab, "sequence_table.txt")'
    shell.write(seq_out)
    shell.write("\n")
    ###############################################################################
    # as we are not using SILVA for this project, this from here doesnt make sense
    # to do.
    # Assign taxonomy
    # tax_out = 'tax <- assignTaxonomy(seqtab, "%s", multithread=TRUE)' % database
    # shell.write(tax_out)
    # shell.write("\n")

    # Write to disk
    # table_out = 'saveRDS(seqtab, "%s")\n' % (os.path.join(working_dir,
                                                          # "dada2",
                                                          # "seqtab_final.rds"))
    # shell.write(table_out)
    # final_table = 'saveRDS(tax, "%s")\n' % (os.path.join(working_dir,
                                                        # "dada2",
                                                        # "tax_final.rds"))
    # shell.write(final_table)
    # run the rscript
    # shell.write("dev.off()\n")
    # shell.write('quit(save = "no", status = 0, runLast = TRUE)\n')

    shell.close()  # MUST do this before try to run it!

    subprocess.check_call([os.path.join("/usr", "bin", "Rscript"),
                           os.path.join(os.getcwd(), "dada2.R")])
    # it might seem mad to trim the seq here, but if doing ITS the seq is too short
    # and the 5.8S region comes through from the left read. So Hvae to trim this off
    # in the final fasta!
    write_as_fasta("sequence_table.txt", right_trim, READ_PREFIX + "_DADA2.fasta")
