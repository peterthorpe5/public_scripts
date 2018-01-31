#!/usr/bin/env python
#
# Tools for working with samtool index
#
# samtools: https://github.com/samtools/samtools
# index allow fast bam file working.
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import subprocess
from collections import namedtuple

from .tools import is_exe, NotExecutableError

# factory class for samtools class returned values
# the order of the outfiles is defined in the build_command self._outfnames

# bam.bai - this is the bam index file
# stderr
Results = namedtuple("Results", "command bai " +
                     "stdout stderr")


class Samtools_IndexError(Exception):
    """Exception raised when Samtools_Index fails"""
    def __init__(self, message):
        self.message = message


class Samtools_Index(object):
    """Class for working with Samtools_index"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, infname):
        """Run Samtools_index on the passed file"""
        self.__build_cmd(infname)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        results = Results(self._cmd, self._outfname, pipe.stdout,
                          pipe.stderr)
        return results

    def __build_cmd(self, infname):
        """Build a command-line for Samtools_index"""
        self._outfname = infname + ".bai"
        cmd = ["samtools",
               "index",
               infname]
        self._cmd = ' '.join(cmd)
