#!/usr/bin/env python
#
# Tools for working with MUSCLE
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import subprocess
from collections import namedtuple

from .tools import is_exe, NotExecutableError

# factory class for Muscle class returned values
# the order of the outfiles is defined in the build_command self._outfnames

# outfile - this is the out file
# stderr
Results = namedtuple("Results",
                     "command outfile stdout stderr")


# TO DO extra options for muscle
class MuscleError(Exception):
    """Exception raised when Muscle fails"""
    def __init__(self, message):
        self.message = message


class Muscle(object):
    """Class for working with MUSCLE"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path


    def run(self, infile, outfile=None, dry_run=False):
        """Run MUSCLE on the single passed file

        Writes the alignment result alongside the input file

        Returns a tuple of output file, and the STDOUT, STDERR returned by the
        muscle run.
        """
        # Construct command and return if a dry run
        self.__build_cmd(infile, outfile)
        if dry_run:
            results = Results(self._cmd, self._outfilename, None, None)
        else:
            pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
            results = Results(self._cmd, self._outfilename, pipe.stdout,
                              pipe.stderr)
        return results

    
    def __build_cmd(self, infile, outfile):
        """Build a command-line for MUSCLE"""
        # Default to a new alignment alongside the input if no output file is
        # specified.
        if outfile is None:
            self._outfilename = os.path.join(os.path.splitext(infile)[0] +
                                             '_muscle.aln')
        else:
            self._outfilename = outfile

        cmd = ["muscle",
               "-in", infile,
               "-out", self._outfilename]
        self._cmd = ' '.join(cmd)
