#!/usr/bin/env python
#
# error_correction * (correct Illumina specific errors)
# http://bioinf.spbau.ru/spades/bayeshammer
# This comes bunduled with SPAdes.
# http://bioinf.spbau.ru/en/content/spades-download-0
# tested with 3.9.0
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import subprocess

from collections import namedtuple

from .tools import is_exe, NotExecutableError

# factory class for EC class returned values
# the order of the outfiles is defined in the build_command self._outfnames

# Left_read_correct - left corrected reads
# right_read_correct - right_read_correct
# unpaired - unpaired.
# stdout
# stderr

Results = namedtuple("Results", "command Left_read_correct " +
                     "right_read_correct unpaired stdout stderr")


class Error_CorrectionError(Exception):
    """Exception raised when flash fails"""
    def __init__(self, message):
        self.message = message


class Error_Correction(object):
    """Class for working with Bayes hammer"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, lreads, rreads, threads, outdir, dry_run=False):
        """Run SPAdes on the read files"""
        assert(lreads != rreads)  # We don't want to EC the same file
        self.__build_cmd(lreads, rreads,
                         threads, outdir)
        if dry_run:
            return(self._cmd)
        # folder checking an making staying in for now
        if not os.path.exists(self._outdirname):
            os.makedirs(self._outdirname)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        results = Results(self._cmd, *self._outfnames, pipe.stdout,
                          pipe.stderr)
        return results

    def __build_cmd(self, lreads, rreads,
                    threads, outdir):
        """Build a command-line for SPAdes to error correct the
         paired end reads.

         Spades has a build in baysian error correction method
         whcih is easy to install.
         it take in
         -l left reads
         -r right_reads
         -threads
         -o out directory.

         the output files name will be names as such: e.g.
         in:  DNAMIX_S95_L001_R1_001.fastq.gz
         out: DNAMIX_S95_L001_R1_001.fastq.00.0_0.cor.fastq.gz
                                     -------------------------

         """
        prefix = os.path.split(lreads)[-1].split("_R")[0]
        L_prefix = os.path.split(lreads)[-1].split(".gz")[0]
        R_prefix = os.path.split(rreads)[-1].split(".gz")[0]
        self._outdirname = os.path.join(outdir)
        self._outfnames = [os.path.join(outdir, "corrected/") +
                           suffix for suffix in
                           (L_prefix + '.00.0_0.cor.fastq.gz',
                           R_prefix + '.00.0_0.cor.fastq.gz',
                           prefix + 'R_unpaired' + '.00.0_0.cor.fastq.gz')]
        cmd = ["spades.py",
               "-1", lreads,
               "-2", rreads,
               "--only-error-correction",
               "--threads", str(threads),
               "-o", self._outdirname]
        self._cmd = ' '.join(cmd)
