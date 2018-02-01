#!/usr/bin/env python
#
# Tools for working with samtools sort
#
# samtools: https://github.com/samtools/samtools
# sorting is required for some tools
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import subprocess
from collections import namedtuple

from .tools import is_exe, NotExecutableError

# factory class for samtools class returned values
# the order of the outfiles is defined in the build_command self._outfnames

# sorted_bam - this is the sorted out bam file.
# stderr
Results = namedtuple("Results", "command sorted_bam stdout stderr")


class Samtools_SortError(Exception):
    """Exception raised when Samtools_Sort fails"""
    def __init__(self, message):
        self.message = message


class Samtools_Sort(object):
    """Class for working with Samtools_Sort"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, infname, outfile, threads):
        """Run Samtools_Sort on the passed file"""
        self.__build_cmd(infname, outfile, threads)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        results = Results(self._cmd, self._outfname, pipe.stdout,
                          pipe.stderr)
        return results

    def __build_cmd(self, infname, outfile, threads):
        """Build a command-line for Samtools_index"""
        self._outfname = outfile
        cmd = ["samtools",
               "sort",
               "-@",
               threads,
               infname,
               '-o', self._outfname]
        self._cmd = ' '.join(cmd)
