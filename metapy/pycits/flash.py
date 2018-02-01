#!/usr/bin/env python
#
# FLASH * (assemble overlapping reads)
# https://ccb.jhu.edu/software/FLASH/
# follow this link to get the binaries.
# https://sourceforge.net/projects/flashpage/
# ?source=typ_redirect
# This has been tested with FLASH-1.2.11

#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

import os
import subprocess

from collections import namedtuple

from .tools import is_exe, NotExecutableError

# factory class for Flash class returned values
# the order of the outfiles is defined in the build_command self._outfnames

# outfileextended - this is the assembled reads
# outfilehist - two column data for assesmbled lengths
# outfilehistogram - ascii histogram
# outfilenotcombined_1 - left reads not used
# outfilenotcombined_2 - right reads not used
# stdout
# stderr
Results = namedtuple("Results", "command outfileextended outfilehist " +
                     "outfilehistogram outfilenotcombined_1 " +
                     "outfilenotcombined_2 stdout stderr")


# To Do pass optional argument to the program.

class FlashError(Exception):
    """Exception raised when flash fails"""
    def __init__(self, message):
        self.message = message


class Flash(object):
    """class for working with paired end read assembly tool Flash"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

    def run(self, lreads, rreads, threads, outdir, prefix, dry_run=False):
        """Run Flash to merge passed read files

        - lreads    - forward reads
        - rreads    - reverse reads
        - threads   - number of threads for flash to use
        - outdir    - output directory for merged output
        - prefix    - file prefix for flash output
        - logger    - stream to write messages
        - dry_run   - if True, returns cmd-line but does not run

        Returns a tuple of output filenames, and the STOUT returned by the
        flash run.
        """
        assert(lreads != rreads)  # We don't want to merge the same file
        self.__build_cmd(lreads, rreads, threads, outdir, prefix)
        if dry_run:
            return(self._cmd)
        # checking and making the output folder is staying in for now
        if not os.path.exists(self._outdirname):
            os.makedirs(self._outdirname)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        results = Results(self._cmd, *self._outfnames, pipe.stdout,
                          pipe.stderr)
        return results

    def __build_cmd(self, lreads, rreads, threads, outdir, prefix):
        """Build a command-line for flash.

        flash takes a path to an output directory PLUS the prefix
        (e.g. temp) of the
        files to write, such that

        -o a/b/temp

        writes files

        a/b/temp.
            temp.extendedFrags.fastq
            temp.hist
            temp.histogram
            temp.notCombined_1.fastq
            temp.notCombined_2.fastq

        """
        self._outfnames = [os.path.join(outdir, prefix) + suffix for suffix in
                           ('.extendedFrags.fastq', '.hist', '.histogram',
                            '.notCombined_1.fastq',
                            '.notCombined_2.fastq')]

        self._outdirname = os.path.join(outdir, "%s_flash" % (prefix))
        cmd = ["flash",
               "--max-overlap", "250",
               "--threads {0}".format(threads),
               "-o", os.path.join(outdir, prefix),
               lreads, rreads]
        self._cmd = ' '.join(cmd)
