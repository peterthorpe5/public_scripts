#!/usr/bin/env python
#
#   module to clean up given files using linux rm command
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import sys
import subprocess
from .tools import is_exe, NotExecutableError
from collections import namedtuple


Results = namedtuple("Results", "command stdout stderr")


class Clean_up(Exception):
    """Exception raised when flash fails"""
    def __init__(self, message):
        self.message = message


class Clean_up(object):
    """Class for cleaning up files."""
    def __init__(self, exe_path):
        """Instantiate with location of executable (rm)"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, infile):
        """Run clean up on unwanted files
        class called linux rm command"""
        self.__build_cmd(infile)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        results = Results(self._cmd, pipe.stdout,
                          pipe.stderr)
        return results

    def __build_cmd(self, infile):
        """Build a command-line to clean up files """
        cmd = ["rm",
               infile]
        self._cmd = ' '.join(cmd)
