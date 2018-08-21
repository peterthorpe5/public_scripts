#!/usr/bin/env python
# title write a shell to run metapy on mutli fastaq files.

import os
import shutil
import subprocess
import sys
 
list_of_read_files = []
seen_set = set([])
# collect all the .gz fastq file names
for filename in os.listdir(".") :
    if not filename.endswith(".gz"):
        continue
    list_of_read_files.append(filename)

# collect all the fastq file names
for filename in os.listdir(".") :
    if not filename.endswith(".fastq"):
        continue
    list_of_read_files.append(filename)

for i in sorted(list_of_read_files):
	prefix = i.split("_R")[0]
	if prefix not in seen_set:
		#print prefix
		seen_set.add(prefix)

# add for different db: "-d", "name.fasta"
for i in sorted(seen_set):
    left = i + "_R1_001.fastq.gz"
    right = i + "_R2_001.fastq.gz"
    commmand = " ".join(["metapy.py",
                         "-l", left,
                         "-r", right,
                         "--thread",
                         "16"])#,"--cleanup", "NO"])
    print(commmand)
    pipe = subprocess.run(commmand, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    print("finished %s" % filename)
    READ_PREFIX = os.path.split(os.path.abspath(left))[-1].split("_R")[0]
    remove_list = ["temp.fasta",
                   "temp.txt",
                   "assembled_fa_and_OTU_db.fasta",
                   "assembled_reads_and_OTU_db.fasta",
                   "db_old_to_new_names.txt",
                   "db_old_to_new_names_vsearch.txt",
                   "assembled_fa_and_OTU_db_vesearch.fasta",
                   "OTU.1.bt2",
                   "OTU.2.bt2",
                   "OTU.3.bt2",
                   "OTU.4.bt2",
                   "OTU.rev.1.bt2",
                   "OTU.rev.2.bt2",
                   "error.log",
                   "dada2_seq_and_OTU_db.fasta",
                   "dada2.R",
                   READ_PREFIX + "_R1_001_primers_trimmed.fastq.gz",
                   READ_PREFIX + "_R2_001_primers_trimmed.fastq.gz",
                   "sequence_table.txt",
                   "temp"]
    for unwanted in remove_list:
        try:
            # pass
            os.remove(unwanted)
        except:
            pass
    try:
        shutil.rmtree(READ_PREFIX)
    except:
        pass
    try:
        shutil.rmtree("dada2")
    except:
        pass
    try:
        shutil.rmtree("filtered")
    except:
        pass


