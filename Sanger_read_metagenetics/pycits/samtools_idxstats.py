#!/usr/bin/env python
#
# Tools for working with samtools idxstats
#
# samtools: https://github.com/samtools/samtools
# idxstats get the number of reads that mapped to
# the database enteries
# - by samtools idxstats.
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import subprocess
from collections import namedtuple

from .tools import is_exe, NotExecutableError

# factory class for samtools class returned values
# the order of the outfiles is defined in the build_command self._outfnames

# cov - this is the reads that map to the database enteries
# - by samtools idxstats.
# stderr
Results = namedtuple("Results", "command cov stdout stderr")


class Samtools_idxstatsError(Exception):
    """Exception raised when samtools idxstats fails"""
    def __init__(self, message):
        self.message = message


class Samtools_Idxstats(object):
    """Class for working with Samtools_idxstats"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, infname, outfile):
        """Run Samtools_Sort on the passed file"""
        self.__build_cmd(infname, outfile)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        results = Results(self._cmd, self._outfname, pipe.stdout,
                          pipe.stderr)
        return results

    def __build_cmd(self, infname, outfile):
        """Build a command-line for Samtools_index"""
        self._outfname = outfile
        # pipe the output through sort and grep without the "8" sign
        cmd = ["samtools", "idxstats", infname,
               "|", "grep", "-v", "'*'",
               "|", "sort", "--reverse", "-n", "-k3",
               ">", self._outfname]
        self._cmd = ' '.join(cmd)
