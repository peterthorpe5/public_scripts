#!/usr/bin/env python
#
# Tools for working with Bowtie2 map
#
# http://bowtie-bio.sourceforge.net/bowtie2/
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import subprocess
from collections import namedtuple
from .tools import is_exe, NotExecutableError

# factory class for bowtie build class returned values
# the order of the outfiles is defined in the build_command self._outfnames

# command -  the command used for the mapping
# sam - this is the outputed sam file.
# stderr
Results = namedtuple("Results", "command sam stdout stderr")


class Bowtie2_MapError(Exception):
    """Exception raised when bowtie2-map fails"""
    def __init__(self, message):
        self.message = message


class Bowtie2_Map(object):
    """Class for working with Bowtie2_Map"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, reads, indexstem, outfilename, threads, fasta=False,
            dry_run=False):
        """Construct and execute a bowtie2 command-line

        reads       - can be a fasta file or a string of left and right reads.
        infnames    - the fasta to index
        indexstem   - the out index prefix
        outdir      - directory for output
        threads     - number of threads to use
        fasta       - Boolean indicating whether reads are FASTA or FASTQ

        if perfect mapping is required consider:
        --score-min 'C,0,-1'
        """
        self.__build_cmd(reads, indexstem, outfilename, threads, fasta)
        print(self._cmd)
        if dry_run:
            results = Results(self._cmd, self._outfname, None, None)
        else:
            pipe = subprocess.run(self._cmd, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True)
            results = Results(self._cmd, self._outfname, pipe.stdout,
                              pipe.stderr)
        return results

    def __build_cmd(self, reads, indexstem, outfilename, threads, fasta):
        """Build a command-line for bowtie2

        For the mapping we will ask for very sensitive alignment, and
        filter .sam output for the number of mismatches we are interested in

        Command-line choices are explained below:

        --very-sensitive    -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
        Some forums say that -N is number of mismatches. THIS IS WRONG
        --no-unal    do not report unaligned reads. Keep the size of the file
        down.
        -p     Threads
        -x     this is the index prefix which must be set up before this is
        run
        -f/-r    -f is a fasta file, -r are comma separated read lists.
        The script currently guesses if it is a fasta ot fq file.
        -S     this is the output sam file. script return the sam file as
        the fast or fq name up to the file extension and the prefix for the
        database searched.
           eg.
           -x targets -f million_pound_gene_db.fasta -S would be ...
           -S million_pound_gene_db_Vs_targets.sam

        We will filter the output on the number of mismatches:
        XM:i:<n> The number of mismatches in the alignment
        """
        self._outfname = outfilename
        if fasta:
            filetype = '-f'
        else:
            filetype = '-q'
        cmd = ["bowtie2",
               "--very-sensitive",
               "--no-unal",
               "-p", threads,
               "-x", indexstem,
               filetype, reads,
               "-S", self._outfname]
        self._cmd = ' '.join(cmd)
