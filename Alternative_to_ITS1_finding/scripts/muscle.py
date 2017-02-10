#!/usr/bin/env python
#
# Tools for working with MUSCLE
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import sys

import subprocess
from collections import namedtuple

from .tools import is_exe, NotExecutableError

# factory class for Muscle class returned values
# the order of the outfiles is defined in the build_command self._outfnames

# outfile - this is the out file
# stderr
Results = namedtuple("Results", "command muscle_outfile " +
                     "stdout stderr")

Results_refine = namedtuple("Results", "command muscle_outfile " +
                     "stdout stderr")


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

    def run(self, infile, outfile, outfolder, dry_run=False):
        """Run MUSCLE on the single passed file
        take the infile as a fasta, the given outfile name and the
        given outfolder name.

        -in infile
        -out outfile
        Returns a tuple of output file, and the STOUT returned by the
        muscle run.
        """

        self._outfolder = os.path.join(outfolder)
        if not os.path.exists(self._outfolder):
            os.makedirs(self._outfolder)

        # ensure the file format is correct - has something weird happened?
        assert (infile.endswith(".fa") or infile.endswith(".fasta"))
        self.__build_cmd(infile, outfile, outfolder)
        if dry_run:
            return(self._cmd)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        results = Results(self._cmd, self._outfilename, pipe.stdout,
                          pipe.stderr)
        return results

    def __build_cmd(self, infile, outfile, outfolder):
        """Build a command-line for MUSCLE"""
        self._outfilename = os.path.join(self._outfolder, outfile)

        cmd = ["muscle",
               "-in", infile,
               "-out",
               self._outfilename]

        self._cmd = ' '.join(cmd)


    def run_refine(self, infile, outfile, outfolder, dry_run=False):
        """Run MUSCLE on the single passed file
        take the infile as a fasta, the given outfile name and the
        given outfolder name.

        -in infile
        -out outfile
        -refine extra_command
        Returns a tuple of output file, and the STOUT returned by the
        muscle run.
        """

        self._outfolder = os.path.join(outfolder)
        if not os.path.exists(self._outfolder):
            os.makedirs(self._outfolder)
        self.__build_cmd_refine(infile, outfile, outfolder)
        if dry_run:
            return(self._cmd)
        # IMPORTANT: check=False here due to problem NULs
        # in files. Bug in muscle?
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=False)
        results = Results(self._cmd, self._outfilename, pipe.stdout,
                          pipe.stderr)
        return results

    def __build_cmd_refine(self, infile, outfile, outfolder):
        """Build a command-line for MUSCLE"""
        self._outfilename = os.path.join(self._outfolder, outfile)

        cmd = ["muscle",
               "-in", infile,
               "-out",
               self._outfilename,
               "-refine"]

        self._cmd = ' '.join(cmd)
