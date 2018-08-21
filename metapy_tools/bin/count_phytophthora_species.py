#!/usr/bin/env python
# title write a shelll to run metapy on mutli fastaq files.

import os
import sys
from collections import defaultdict

SPECIES_DIC_COUNT = defaultdict(int)
SPECIES_READ_COUNT = defaultdict(int)

FILECOUNT = 0

def test_line(line):
    """returns true lines. Not comments or blank line"""
    if not line.strip():
        return False  # if the last line is blank
    if line.startswith("#"):
        return False  # comment line
    if line.startswith("    #"):
        return False  # comment line
    return line


for filename in os.listdir("."):
    if not filename.endswith(".RESULTS"):
        continue
    with open(filename) as handle:
        FILECOUNT += 1
        for line in handle:
            if test_line(line):
                cluster, species, reads = line.split("\t")
                SPECIES_DIC_COUNT[species] += 1
                SPECIES_READ_COUNT[species] += int(reads)

print("\n%d files interrogated" % FILECOUNT)

print("\n\tthe following Phytophthoras/controls were found")
print(SPECIES_DIC_COUNT)

print("\nif we divide those identified  by the number of files:")

for keys, vals in SPECIES_DIC_COUNT.items():
    print("%s\t%.2f percent" % (keys, 100*(float(vals)/FILECOUNT)))

print("\n\tof these - how many reads hit the species over all the file?")
print(SPECIES_READ_COUNT)

print("\n\tif we divide those reads by the number of files:")

for keys, vals in SPECIES_READ_COUNT.items():
    print("%s\t%.2f\treads per file" % (keys, (float(vals)/FILECOUNT)))

                
