#!/usr/bin/env python
#
# Tools for working with seq_crumbs scripts
#
# seq_crumbs: https://github.com/JoseBlanca/seq_crumbs
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard

import os
import sys

from subprocess import call
from .tools import is_exe


class Trim_Quality(object):
    """Class for working with trim_quality"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("No trim_quality script available (exiting)")
            sys.exit(1)
        self._exe_path = exe_path
        self.format = 'fastq'

    def run(self, infname, outdir):
        """Run trim_quality on the passed file"""
        self.__build_cmd(infname, outdir)
        msg = ["Running...", "\t%s" % self._cmd]
        for m in msg:
            self._logger.info(m)
        retcode = call(self._cmd, shell=True)
        if retcode < 0:
            self._logger.error("trim_quality terminated by " +
                               "signal %s" % -retcode)
            sys.exit(1)
        else:
            self._logger.info("trim_quality returned %s" % retcode)
        return self._outfname

    def __build_cmd(self, infname, outdir):
        """Build a command-line for trim_quality"""
        fname = os.path.split(infname)[-1]
        outfname = fname.split('.fastq')[0] + '_trimmed' + '.' + self.format
        outfilename = os.path.join(outdir, outfname)
        cmd = ["trim_quality",
               "-t", self.format,
               "--paired_reads", infname,
               "-o", outfilename]
        self._cmd = ' '.join(cmd)
        self._outfname = outfilename


class Convert_Format(object):
    """Class for working with convert_format"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("No convert_format script available (exiting)")
            sys.exit(1)
        self._exe_path = exe_path
        self.informat = "fastq"
        self.outformat = "fasta"

    def run(self, infname, outdir, logger=None):
        """Run convert_format on the passed file"""
        self._outdirname = os.path.join(outdir)
        if not os.path.exists(self._outdirname):
            if logger:
                self._logger.info("Creating output directory: %s" %
                                  self._outdirname)
                os.makedirs(self._outdirname)
        self.__build_cmd(infname, self._outdirname)
        msg = ["Running...", "\t%s" % self._cmd]
        for m in msg:
            self._logger.info(m)
        retcode = call(self._cmd, shell=True)
        if retcode < 0:
            self._logger.error("convert_format terminated by " +
                               "signal %s" % -retcode)
            sys.exit(1)
        else:
            self._logger.info("convert_format returned %s" % retcode)
        return self._outfname

    def __build_cmd(self, infname, outdir):
        """Build a command-line for convert_format"""
        fname = os.path.split(infname)[-1]
        outfname = os.path.splitext(fname)[0] + ".%s" % self.outformat
        outfilename = os.path.join(outdir, outfname)
        cmd = ["convert_format",
               "-t", self.informat,
               "-o", outfilename,
               "-f", self.outformat,               
               infname]
        self._cmd = ' '.join(cmd)
        self._outfname = outfilename
