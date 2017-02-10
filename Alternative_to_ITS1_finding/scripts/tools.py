#!/usr/bin/env python
#
#   module to clean up given files using linux rm command
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import sys
from subprocess import check_output, CalledProcessError


class NotExecutableError(Exception):
    """Exception raised when expected executable is not executable"""
    def __init__(self, message):
        self.message = message


def is_exe(filename):
    """Returns True if path is to an executable file"""
    if os.path.isfile(filename) and os.access(filename, os.X_OK):
        return True
    else:
        try:
            exefile = check_output(["which", filename]).strip()
        except CalledProcessError:
            raise NotExecutableError("{0} does not exist".format(filename))
    return os.path.isfile(exefile) and os.access(exefile, os.X_OK)
