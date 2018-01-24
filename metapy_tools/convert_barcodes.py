#!/usr/bin/env python
# title write a shelll to run metapy on mutli fastaq files.

import os
import sys
from sys import stdin,argv
from Bio import SeqIO

wanted = set("""
M01157:20:000000000-BGPHF:1:1102:19601:2002	M01157:20:000000000-BGPHF:1:1102:19602:2026	M01157:20:000000000-BGPHF:1:1112:8977:15172	M01157:20:000000000-BGPHF:1:2114:5653:9995	M01157:20:000000000-BGPHF:1:2111:21246:2863	M01157:20:000000000-BGPHF:1:2112:10983:20935	M01157:20:000000000-BGPHF:1:1102:6467:13977	M01157:20:000000000-BGPHF:1:1109:13351:20528	M01157:20:000000000-BGPHF:1:1109:13855:16149	

""".split())

for record in SeqIO.parse(argv[1], "fastq"):
    if record.id in wanted:
        print("%s\t%s" % (record.id, record.seq))
        print(record.letter_annotations["phred_quality"])
