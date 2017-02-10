#!/usr/bin/env python
#
# Tools for working with Gblocks
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import sys

import subprocess
from collections import namedtuple

from .tools import is_exe, NotExecutableError

# factory class for Gblocks class returned values
# the order of the outfiles is defined in the build_command self._outfnames

# outfile - this is the out file
# stderr
Results = namedtuple("Results", "command Gblocks_outfile " +
                     "stdout stderr")

Results_refine = namedtuple("Results", "command Gblocks_outfile " +
                     "stdout stderr")


# TO DO extra options for Gblocks
class GblocksError(Exception):
    """Exception raised when Gblocks fails"""
    def __init__(self, message):
        self.message = message


class Gblocks(object):
    """Class for working with Gblocks"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, infile, outprefix, dry_run=False):
        """Run Gblocks on the single passed file
        take the infile as a fasta, the outfile will
        be infile + prefix
        http://molevol.cmima.csic.es/castresana/Gblocks/
        Gblocks_documentation.html#command_line
        """
        # ensure the file format is correct
        # - has something weird happened?
        assert (infile.endswith(".fa") or infile.endswith(".fasta"))
        self.__build_cmd(infile, outprefix)
        if dry_run:
            return(self._cmd)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        results = Results(self._cmd, self._outfilename, pipe.stdout,
                          pipe.stderr)
        return results

    def __build_cmd(self, infile, outprefix):
        """Build a command-line for Gblocks"""
        self._outfilename = os.path.join(infile + outprefix)
        cmd = ["Gblocks",
               infile,
               "-t=d",
               "-b3=200",
               "-b4=200",
               "-b5=a",
               "-e="
               outprefix]]

        self._cmd = ' '.join(cmd)
