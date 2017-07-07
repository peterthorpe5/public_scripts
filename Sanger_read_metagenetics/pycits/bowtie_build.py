#!/usr/bin/env python
#
# Tools for working with Bowtie2 build
#
# http://bowtie-bio.sourceforge.net/bowtie2/
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import subprocess
from collections import namedtuple

from .tools import is_exe, NotExecutableError

# factory class for bowtie build class returned values
# the order of the outfiles is defined in the build_command self._outfnames

# index - this is the index file generated.
# stderr
Results = namedtuple("Results", "command index stdout stderr")


class Bowtie2_BuildError(Exception):
    """Exception raised when bowtie2-build fails"""
    def __init__(self, message):
        self.message = message


class Bowtie2_Build(object):
    """Class for working with bowtie2-build"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, infname, outstem, dry_run=False):
        """Construct and execute a bowtie2-build command-line"""
        self.__build_cmd(infname, outstem)
        if dry_run:
            results = Results(self._cmd, self._outfname, None, None)
        else:
            pipe = subprocess.run(self._cmd, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True)
            results = Results(self._cmd, self._outfname, pipe.stdout,
                              pipe.stderr)
        return results

    def __build_cmd(self, infname, outstem):
        """Build a command-line for bowtie2-build"""
        self._outfname = outstem
        cmd = ["bowtie2-build",
               "--quiet",
               "-f",
               infname,
               self._outfname]
        self._cmd = ' '.join(cmd)
