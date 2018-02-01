#!/usr/bin/env python
#
# Tools for working with FastQC
#
# FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import subprocess
from collections import namedtuple

from .tools import is_exe, NotExecutableError

# factory class for Flash class returned values
# the order of the outfiles is defined in the build_command self._outfnames

# fastqc_html - this is the html file
# fastqc_zip -  this is the output zip folder.
# stderr
Results = namedtuple("Results",
                     "command htmlfile zipfile stdout stderr")


class FastQCError(Exception):
    """Exception raised when fastqc fails"""
    def __init__(self, message):
        self.message = message


class FastQC(object):
    """Class for working with FastQC"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, infnames, outdir, dry_run=False):
        """Run fastqc on the passed file"""
        self.__build_cmd(infnames, outdir)
        if not os.path.exists(self._outdirname):
            os.makedirs(self._outdirname)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        results = Results(self._cmd, *self._outfnames, pipe.stdout,
                          pipe.stderr)
        return results

    def __build_cmd(self, infname, outdir):
        """Build a command-line for fastqc
        fastqc returns a html and a .zip folder after the given read name
        """
        prefix = os.path.split(infname)[-1].split(".fa")[0]
        self._outfnames = [os.path.join(outdir, prefix) + suffix for suffix in
                           ('_fastqc.html', '_fastqc.zip')]
        self._outdirname = os.path.join(outdir)
        cmd = ["fastqc",
               infname,
               "-o", self._outdirname]
        self._cmd = ' '.join(cmd)
