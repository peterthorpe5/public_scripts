#!/usr/bin/env python
#
# cd_hit * (clusterassembled reads with database)
# cd_hit_est used as this is the nt clustering tool
# http://weizhongli-lab.org/lab-wiki/doku.php?id=cd-hit-user-guide
# follow this link to get the download.
# https://github.com/weizhongli/cdhit
# cd_hit-0.9.10-bin-64.tar.gz
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import subprocess

from collections import namedtuple

from .tools import is_exe, NotExecutableError


# factory class for cd_hit class returned values
Results = namedtuple("Results", "command fastaout clusters " +
                     "stdout stderr")


class Cd_hit_Error(Exception):
    """Exception raised when cd_hit fails"""
    def __init__(self, message):
        self.message = message


class Cd_hit(object):
    """Class for working with cd_hit"""

    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, fasta_in, threads, threshold, outdir, prefix, dry_run=False):
        """Run cd_hit to cluster passed fasta files

        - fasta_in    - fasta file to be clustered
        - threshold - threshold to cluster at
        - threads   - number of threads for cd_hit to use
        - outdir    - output directory for clustering  output
        - prefix    - file prefix for cd_hit output
        - dry_run   - if True, returns cmd-line but does not run

        Returns a tuple of output filenames, and the STOUT returned by the
        cd_hit run.
        """
        self.__build_cmd(fasta_in, threads, threshold, outdir, prefix)
        if dry_run:
            return(self._cmd)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)

        results = Results(self._cmd, *self._outfnames, pipe.stdout,
                          pipe.stderr)
        return results

    def __build_cmd(self, fasta_in, threads, threshold, outdir, prefix):
        """Build a command-line for cd_hit_est.

        cd_hit takes a path to an output directory PLUS the prefix of the
        files to write, such that

        -o a/b/cdefg

        writes files
        a/b/cdefg
        a/b/cdefg.clstr

        and so on.

        -d added to the command is so the output clusters will write out the
        names up to 500 letters long. The default chops these at 20.
        (too short)

        -M added to allow unlimited memeroy - not a problem for
        small jobs. If job are big, we will have to alter this.
        """
        # outfiles are name WhatEver.out + .bak.clstr and + .clstr
        self._outfnames = [os.path.join(outdir, prefix) + suffix for suffix in
                           ('.fasta', '.clstr')]
        cmd = ["cd-hit-est",
               "-i", fasta_in,
               "-o", os.path.join(outdir, prefix),
               "-T {0}".format(threads),
               "-M", "0",
               "-c", str(threshold),
               "-d", "500"]
        self._cmd = ' '.join(cmd)
