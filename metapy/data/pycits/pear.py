#!/usr/bin/env python
#
# PEAR * (assemble overlapping reads)
# https://github.com/xflouris/PEAR
# follow this link to get the binaries.
# http://sco.h-its.org/exelixis/web/software/pear/files/
# pear-0.9.10-bin-64.tar.gz
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

# WARNING: this module is not yet tested AT ALL.

import os
import subprocess

from collections import namedtuple

from .tools import is_exe, NotExecutableError


# factory class for Pear class returned values
Results = namedtuple("Results", "command outfileassembled outfilediscarded " +
                     "outfileunassmbledfwd outfileunassembledrev " +
                     "stdout stderr")


class PearError(Exception):
    """Exception raised when pear fails"""
    def __init__(self, message):
        self.message = message


class Pear(object):
    """Class for working with PEAR"""

    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, lreads, rreads, threads, outdir, prefix, dry_run=False):
        """Run PEAR to merge passed read files

        - lreads    - forward reads
        - rreads    - reverse reads
        - threads   - number of threads for pear to use
        - outdir    - output directory for merged output
        - prefix    - file prefix for pear output
        - logger    - stream to write messages
        - dry_run   - if True, returns cmd-line but does not run

        Returns a tuple of output filenames, and the STOUT returned by the
        pear run.
        """
        assert(lreads != rreads)  # We don't want to merge the same file
        self.__build_cmd(lreads, rreads, threads, outdir, prefix)
        if dry_run:
            return(self._cmd)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        results = Results(self._cmd, *self._outfnames, pipe.stdout,
                          pipe.stderr)
        return results

    def __build_cmd(self, lreads, rreads, threads, outdir, prefix):
        """Build a command-line for pear.

        pear takes a path to an output directory PLUS the prefix of the
        files to write, such that

        -o a/b/cdefg

        writes files

        a/b/cdefg.assembled.fastq
        a/b/cdefg.discarded.fastq

        and so on.
        """
        self._outfnames = [os.path.join(outdir, prefix) + suffix for suffix in
                           ('.assembled.fastq', '.discarded.fastq',
                            '.unassembled.forward.fastq',
                            '.unassembled.reverse.fastq')]
        cmd = ["pear",
               "-f", lreads,
               "-r", rreads,
               "--threads {0}".format(threads),
               "-o", os.path.join(outdir, prefix)]
        self._cmd = ' '.join(cmd)
