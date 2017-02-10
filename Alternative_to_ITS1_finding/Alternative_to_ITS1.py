#!/usr/bin/env python
# author: Peter Thorpe and Leighton Pritchard. The James Hutton Insitute,Dundee,UK.
# 2017 Feb
# Title:
# this is a master script to identify an alternative to ITS1 for identification
# of a species.
# BASED on the output Busco Version 1.1b
##############################################################################


import sys
# check this is running under python 3
if sys.version_info[0] < 3:
    # not sure about the error typ ehere?
    raise ImportError ("\n\t\tMust be using Python 3\n")
    sys.exit(1)

import argparse
import os
import warnings
from collections import defaultdict
import subprocess
import logging
import logging.handlers
import traceback
import subprocess
from argparse import ArgumentParser
import time
import gzip


from Bio import SeqIO
from Bio import BiopythonParserWarning

from scripts import muscle
from scripts.convert_GBK_fa import convert_file
from scripts.identify_single_copy_busco_full_table import parse_tab_outfile
from scripts.identify_single_copy_busco_full_table import complete_single_busco
from scripts.write_fasta import write_fasta



HOW = """\n    1) Run Busco on the genomes available.
http://busco.ezlab.org/

    2) parse the full_table_GENOME_NAME_BUSCO file.
 - identify those that are single copy and complete.

    3) Convert the gbk files to fasta in the gb folder,
for those single copy BUSCOS of interest. - exons only??
- currently NOT working.

    3b) get the genomic regions (with Introns) as defined
in the full table.

    4) Identify common single copy busco in above
threshold number of the genomes in question.

    5) alignement: Muscle and refine the muscle alignment
must be in your PATH as muscle,
http://www.drive5.com/muscle/downloads.htm .

    6) run GBLOCKs to filter the alignment,
is this a good idea?:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4538881/
I expect these alignments to be poor due to the size.
Gblocks_0.91b musct be in your PATH as: Gblocks. 
http://molevol.cmima.csic.es/castresana/Gblocks/
Gblocks_documentation.html
ALTER PARAMETERS?
TODO: write a class for this to take parameters. 

    7) convert file to phylip. Biopython

    8) DN/DS ratio: CodonPhyml. 
codonphyml must be in your PATH.
https://sourceforge.net/projects/codonphyml/files/.
https://academic.oup.com/mbe/article/30/6/1270/1134508/.
CodonPhyML-Fast-Maximum-Likelihood-Phylogeny.
http://gyra.ualg.pt/_static/codonPhyML_Manual.pdf.

    9) can primers be disgned for all?.
Having decided how to do this yet.

Requires:
    Biopython,
    muscle,
    Gblock,
    codonphyml (and its dependencies)
"""

###############################################################################

def get_args():
    parser = argparse.ArgumentParser(
    description="Pipeline to try to identify alternative for metabarcoding " +
    ". An Alternative to ITS \n" + HOW,
    add_help=False)
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-t", "--threshold", dest='threshold',
                          action="store", default=0.7,
                          type=float,
                          help="the minimum amount of genomes that BUSCO" +
                          " identified gene must be in, to be considered "
                          "for analysis " +
                          "ranges 0.0 to 1.0")
    optional.add_argument("--min_len", dest='min_len',
                          action="store", default="250",
                          type=str,
                          help="the minimum len Gblock will look for " +
                          "when filtering the alignment. "
                          "The can be altered to suit your seq technology " +
                          "approach")
    optional.add_argument("--AllowedGap", dest='AllowedGap',
                          action="store", default="a",
                          type=str,
                          choices=["a", "h", "n"],
                          help="Allowed Gap positions n (none) " +
                          ", h (half), a (all) . Default a (all)")
    optional.add_argument("--max_non_contig", dest='max_non_contig',
                          action="store", default="8",
                          type=str,
                          help="the max number of non contig position " +
                          "in the alignment when filtering the alignment. ")
    optional.add_argument("--use_genbank", dest='use_genbank',
                          action="store",
                          choices=["True", "False"],
                          help="to get the seq from the gb files " +
                          "outputted by BUSCO. These are often missing " +
                          "default: False", default=False,
                          type=str)
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    args = parser.parse_args()
    return args

###############################################################################

if __name__ == "__main__":
    args = get_args()
    # Set up logging
    logger = logging.getLogger('Alternative_to_ITS1.py: %s' % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)

    try:
        logstream = open("Alternative_to_ITS1_testing.txt", 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        logger.error("Could not open %s for logging" %
                     args.logfile)
        sys.exit(1)

    # Report input arguments
    logger.info("Command-line: %s" % ' '.join(sys.argv))
    logger.info("Starting testing: %s" % time.asctime())
    logger.info(" HOW: %s" % HOW)

    ###########################################################################
    # all the genomes e.g. 0.7 means in 70% of all genomes.
    threshold = args.threshold  # threshold for the num of busco GOI present in
    use_genbank = args.use_genbank  # get the seq from gbk file.
    #Most of these are missing
    ###########################################################################

    # Make folder to put files in

    OUTFOLDER_ALI = "muscle_aligned"
    OUTFOLDER_ALI_REFINE = "muscle_aligned_refined"
    EOG_FASTA_FILE = "EOG_FASTA_FILE"

    folders = [OUTFOLDER_ALI, OUTFOLDER_ALI_REFINE, EOG_FASTA_FILE]
    for i in folders:
        if not os.path.exists(i):
            logger.info("making folder %s" % i)
            os.makedirs(i)


    ###########################################################################

    # 1) script assumes BUSCO has already been run.

    # can generate all the .sh scripts needed for running BUSCO using:
    # ./scripts/write_busco.py
    logger.info("1) script assumes BUSCO has already been run.")

    # create dicts: genome name to prefix used
    genome_prefix = dict()
    # genome.fa = genome
    genome_folder = dict()
    # genome.fa = run_BUSCO_prefix
    genome_full_table = dict()
    Busco_GOI = defaultdict(int)  # busco genes of interests
    genom_scaff_start_stop = defaultdict(list)
    genome_count = 0

    ###########################################################################
    # 2) parse the full tables. to identify busco single copy genes

    logger.info("2) parse the full tables.")

    logger.info("looking for compressed genome seq fna.gz")
    for filename in os.listdir(".") :
        if filename.endswith(".fna.gz"):
            # decompress genomic files
            cmd = " ".join(["gunzip", filename])
            logger.info("decompressing %s" % filename)
            pipe = subprocess.run(cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)

    for filename in os.listdir(".") :
        if not filename.endswith(".fna"):
            continue
        genome_count = genome_count + 1
        logger.info("genome: %s, count %d" % (filename,
                                              genome_count))
        name_prefix = filename.split("_genomic")[0]
        full_table = "full_table_" + name_prefix
        folder = "run_busco_" + name_prefix
        genome_prefix[filename] = name_prefix
        genome_folder[filename] = folder
        genome_full_table[filename] = full_table

    logger.info("there are %d genomes in this folder" %
                genome_count)


    for genome, prefix in genome_prefix.items():
        filepath = os.path.join(genome_folder[genome],
                                genome_full_table[genome] + "_BUSCO")
        # check this file exits. Its placement is weird.
        if not os.path.isfile(filepath):
            # sometime the output is current directory
            filepath = os.path.join(genome_full_table[genome] + "_BUSCO")
        if not os.path.isfile(filepath):
            logger.error("Could'nt open %s busco table for genome %s" %
                         (filepath, genome))
            sys.exit(1)

        # create a genome named based dictionary to assigne scaff start stop to
        # function will populate this
        genom_scaff_start_stop ["%s_scaff_start_stop" % (genome)] = []
        complete_busco, name_list,\
                        genom_scaff_start_stop = \
                        complete_single_busco(filepath,
                                             genome,
                                             genom_scaff_start_stop)
        logger.info("genome %s has %d complete BUSCO EOG" % (genome,
                                                             len(name_list)))
        # iterate through the return names and count them
        for busco in name_list:
            Busco_GOI[busco] += 1

    # to identify those that are in most of the geneoms interogated
    logger.info("only accept EOG in %.1f percent of all genomes" % (
                100*threshold))
    GOI_passed_threshold = []
    for key, vals in Busco_GOI.items():
        ratio_GOI_present = float(vals)/float(genome_count)
        # threshold defines at top of the scripts.
        if ratio_GOI_present >= threshold:
            GOI_passed_threshold.append(key)
    logger.info("%d busco GOI passed threshold" %(len(GOI_passed_threshold)))


    ###########################################################################
    # 3) convert the gbk of the buscos of interest to fa
    # but rename the id to the prefix so we can keep track of them later

    # set to false by default as too many re missing.
    if use_genbank == True:
        logger.info("genbank=true. Careful many of these are missing")
        # interate through the list of passed BUSCOs names
        for GOI in GOI_passed_threshold:
            for genome, prefix in genome_prefix.items():
                gb_name = os.path.join(genome_folder[genome], "gb",
                                       GOI + ".raw.gb")
                fa_out = os.path.join(genome_folder[genome], "gb",
                                              GOI + ".nt.fasta")
                if not os.path.isfile(gb_name):
                    print("WARNING - missing GB file %s" % gb_name)
                    continue
                # convert gbk to fasta AND add prefix to record name
                try:
                    # We expect lots of these warnings:
                    # BiopythonParserWarning: Malformed LOCUS line found -
                    # is this correct?
                    # Using the with statement to mk this filter only temporary
                    logger.info("biopython gbk warning OFF")
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", BiopythonParserWarning)
                        # expect one and only one record!
                        record = SeqIO.read(gb_name, "gb")
                except ValueError:
                    logger.info("ERROR - problem with GB file %s" % gb_name)
                    continue
                record.id = prefix + "_" + GOI.replace("BUSCOa", "")
                record.name == record.id  # to avoid oddity in FASTA output
                count = SeqIO.write(record, fa_out, "fasta")
                assert count == 1, logger.info("GBK has more than one entry")

    ###########################################################################
    # 3b) get the regions of intersts based on genomic DNA, the gb may have exons
    # only. genom_scaff_start_stop has been populated in section 2.
    # the genome is a variable that altered with the genome name.

    # save the sequence to a dictionary so we can get them later
    sequences_dict = defaultdict(list)

    for genome, prefix in genome_prefix.items():
        logger.info("3b) indexing genome %s" % genome)
        genome_seq = SeqIO.index(genome, "fasta")
        coordinates_list = genom_scaff_start_stop [genome]
        for i in coordinates_list:
            # TODO - Check start coord +/- one?
            # 0 computer counting
            BUSCO_group, scaff, start, stop = i.split("\t")
            # break if the busco is not in the passed list
            if not BUSCO_group in GOI_passed_threshold:
                continue
            seq_record = genome_seq[scaff]
            if int(start)> int(stop):
                # need to reverse complement
                seq_record.seq = seq_record.seq[int(start): int(stop)]
                logger.info("start > stop. Rev_complement % s" %
                            BUSCO_group)
                seq_record.seq = seq_record.seq.reverse_complement()
            else:
                seq_record.seq = seq_record.seq[int(start): int(stop)]
            file_name = BUSCO_group + "_genomic_coord.fasta"
            seq_record.id = prefix + "_" + BUSCO_group
            seq_record.description = ""
            #dataformatted = "%s\t%s" % (seq_record.id, str(seq_record.seq))
            sequences_dict[file_name].append(seq_record)
            continue # this works for now
            #SeqIO.write(seq_record, file_name, "fasta")

    # iterate over the dictionay and write to seq to files.
    for filename, records in sequences_dict.items():
        filename = open(filename, "w")
        #print (filename, records)
        for seq_record in records:
            # logger.info("writing to %s %s\n " % (filename, seq_record.id))
            seq_record.description = ""
            SeqIO.write(seq_record, filename, "fasta")
        filename.close()


    ###########################################################################
    # 4b) need to make a single file for each PASSED busco for alignmnet

    # set to false by default as too many re missing.
    if use_genbank == True:
        for GOI in GOI_passed_threshold:
            # interate through the list of passed BUSCOs names
            # open a file to put all the corresponding busco seq to
            GOI_file = open(GOI + ".fasta", "w")

            for genome, prefix in genome_prefix.items():
                fasta = os.path.join(genome_folder[genome], "gb",
                                     GOI + ".nt.fasta")
                if os.path.isfile(fasta):
                    # print ("looking for %s " % fasta)
                    write_fasta(fasta,GOI_file)
                else:
                    logger.info("missing %s" % fasta)
                    continue
                    # possibly a file that was never made GBK files are missing.
            GOI_file.close()
    else:
        logger.info("4b) not using gbk files")

    ###########################################################################
    # 5) align these files

    try:
        # check that muscle is in PATH
        obj = muscle.Muscle("muscle")
    except ValueError:
        logger.info("you must have uscle in your PATH as muscle")
        os._exit()
    logger.info("5) aligning, moving and refining alinments")


    # counter for reporting/ or not!
    report_alin = 0
    report_alin_refine = 0

    # TODO: pass out each alignment to a different processor?

    for filename in os.listdir(".") :
        if not filename.endswith(".fasta"):
            continue
        if "BUSCOaEOG7BKR80" in filename:
            logger.info("%s is being missed due to errors" % filename)
            # reproducibly cause muscle to
            #*** glibc detected *** muscle: double
            # free or corruption (!prev): 0x000000000251f130 ***
            continue
        outfile = filename.split(".fasta")[0] + "_align.fa"
        # run the muscle class
        align = obj.run(filename, outfile, OUTFOLDER_ALI)

        if report_alin == 0:
            logger.info("only reporting 1 command to reduce file size")
            logger.info("align command = %s", align.command)
            report_alin = report_alin + 1
        # move the fast file to fasta folder.
        # more elegent to directly make them there instead
        cmd = " ".join(["mv", filename, "EOG_FASTA_FILE"])
        pipe = subprocess.run(cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        # refine the aligmnet
        outfile_re = filename.split(".fasta")[0] + "_ali_refi.fasta"
        # errors have been found in the muscle aligned file. lines on NUL
        # therefore wran abou these files and move on.
        try:
            align_refine = obj.run_refine(align.muscle_outfile, outfile_re,
                                          OUTFOLDER_ALI_REFINE)
        except:
            # ValueError:  # This is bad form but muscle keeps crashing
            logger.info("WARNING this has been BLINDLY passed!!")
            # IMPORTANT: check=False here due to problem NULs
            # in muscle run_refine class
            logger.info("ERROR - problem with aliment file %s" %
                        align.muscle_outfile)
            # there has been a problems with alignmnet
            logger.info("skipping this file")
            pass
        if report_alin_refine == 0:
            logger.info("only reporting 1 refine command to reduce file size")
            report_alin_refine = report_alin_refine + 1
            logger.info("refine command = %s", align_refine.command)

        # delete the temp aligned file
        cmd = " ".join(["rm", os.path.join(OUTFOLDER_ALI, outfile)])
        pipe = subprocess.run(cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        
    logger.info("5) Alinmnet done")


    ###########################################################################
    # 6)  run GBLOCK
    logger.info("6)run GBLOCK to filter the alignment")

    # filter the alignmnet before phy conversion
    # Gblocks ${f} -t=d -b3=200 -b4=200 -b5=a -e=_gblock.out
    logger.info("Gblocks comands options:")
    logger.info("max_non_contig  = %s " % args.max_non_contig)
    logger.info("min_len = %s " % args.min_len)
    logger.info("AllowedGap = %s " % args.AllowedGap)
    gbloc_count = 0

    # TODO: write a class for this to take parameters
    for filename in os.listdir(OUTFOLDER_ALI_REFINE):
        if not filename.endswith(".fasta"):
            continue
        filename = os.path.join(OUTFOLDER_ALI_REFINE, filename)
        gblock_cmd = " ".join(["Gblocks",
                               filename,
                               "-t=d",
                               "-b3=%s" % args.max_non_contig,
                               "-b4=%s" % args.min_len,
                               "-b5=%s" % args.AllowedGap])
        pipe = subprocess.run(gblock_cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        if gbloc_count == 0:
            logger.info("only reporting one command")
            logger.info("Gblocks command %s " % gblock_cmd)
            gbloc_count = gbloc_count + 1
                
    logger.info("6)Gblocks done")

    ###########################################################################
    # 7)  convert file to phylip

    phylip_count = 0
    logger.info("7) need to convert gblock out to .phy")
    Phylip_files = os.path.join("Phylip_files")

    if not os.path.exists(Phylip_files):
        logger.info("making folder %s" % Phylip_files)
        os.makedirs(Phylip_files)

    for filename in os.listdir("."):
        if not filename.endswith("gb"):
            continue
        alignment = AlignIO.read(open(filename),
                                 "fasta",
                                 alphabet=Gapped(IUPAC.ambiguous_dna))
        outfile = os.path.join(Phylip_files,
                               "%s.phy" %(filename.split("._gblock.out")[0]))
        AlignIO.write(alignment, outfile, "phylip-relaxed")
        if phylip_count == 0:
            logger.info("only reporting one command")
            logger.info("phylip converted %s to %s command %s " % (filename,
                                                                   outfile))
            phylip_count = phylip_count + 1

    # gblocks to pull out regions of interest.
    # filter alignments: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4538881/

    ###########################################################################
    # 7)  run GBLOCK, convert file to phylip
    logger.info("8) DN/DS - codonphyml")
    # codonphyml
    # codonphyml -i ${f} -m GY --fmodel F3X4 -t e -f empirical -w g -a e
    # look up parameters
    for filename in os.listdir(Phylip_files):
        if not filename.endswith(".phy"):
            continue
        # page 116 117 of manual:
        # http://gyra.ualg.pt/_static/codonPhyML_Manual.pdf
        codonphy_cmd = " ".join(["codonphyml",
                                 "-i",
                                 filename,
                                 "-m",
                                 "GY",
                                 "--fmodel",
                                 "F3X4",
                                 "-t",
                                 "e",
                                 "-f",
                                 "empirical",
                                 "-w",
                                 "g",
                                 "-a",
                                 "e"])
        pipe = subprocess.run(codonphy_cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
    logger.info("pipeline up to primer design finished")
    logger.info("compress the genomes to save space")

    for filename in os.listdir(".") :
        if not filename.endswith(".fna"):
            continue
        # decompress genomic files
        cmd = " ".join(["gzip", filename])
        logger.info("compressing %s" % filename)
        pipe = subprocess.run(cmd, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)

