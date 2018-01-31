#!/usr/bin/env python
#
# Tools for working with QIIME scripts
#
# QIIME: http://qiime.org/
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard

import os
import sys

from subprocess import call
from .tools import is_exe


class Pick_Otus(object):
    """Class for working with pick_otus.py"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("No pick_otus.py script (exiting)")
            sys.exit(1)
        self._exe_path = exe_path

    def run(self, infnames, refname, outdir):
        """Run pick_otus.py on the passed file"""
        self.__build_cmd(infnames, refname, outdir)
        msg = ["Running...", "\t%s" % self._cmd]
        for m in msg:
            self._logger.info(m)
        retcode = call(self._cmd, shell=True)
        if retcode < 0:
            self._logger.error("pick_otus.py terminated by " +
                               "signal %s" % -retcode)
            sys.exit(1)
        else:
            self._logger.info("pick_otus.py returned %s" % retcode)
        return self._outdirname

    def __build_cmd(self, infname, refname, outdir):
        """Build a command-line for pick_otus.py"""
        self._outdirname = os.path.join(outdir, "qiime_uclust_OTUs")
        cmd = ["pick_otus.py",
               "-m", "uclust_ref",
               "-s", "0.99",
               "-i", infname,
               "-r", refname,
               "-o", self._outdirname]
        self._cmd = ' '.join(cmd)


class Pick_Closed_Ref_Otus(object):
    """Class for working with pick_closed_reference_otus.py"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("No pick_closed_reference_otus.py script " +
                               "(exiting)")
            sys.exit(1)
        self._exe_path = exe_path

    def run(self, infnames, refname, outdir):
        """Run pick_closed_reference_otus.py on the passed file"""
        self.__build_cmd(infnames, refname, outdir)
        msg = ["Running...", "\t%s" % self._cmd]
        for m in msg:
            self._logger.info(m)
        retcode = call(self._cmd, shell=True)
        if retcode < 0:
            self._logger.error("pick_closed_reference_otus.py terminated by " +
                               "signal %s" % -retcode)
            sys.exit(1)
        else:
            self._logger.info("pick_closed_reference_otus.py returned " +
                              "%s" % retcode)
        return self._outdirname

    def __build_cmd(self, infname, refname, outdir):
        """Build a command-line for pick_closed_reference_otus.py"""
        self._outdirname = os.path.join(outdir, "qiime_closedref_OTUs")
        cmd = ["pick_closed_reference_otus.py", "-f",
               "-i", infname,
               "-r", refname,
               "-o", self._outdirname]
        self._cmd = ' '.join(cmd)


class Join_Paired_Ends(object):
    """Class for working with join_paired_ends.py"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("No join_paired_ends.py script (exiting)")
            sys.exit(1)
        self._exe_path = exe_path

    def run(self, infnames, outdir):
        """Run joined_paired_ends.py on the passed file"""
        self.__build_cmd(infnames, outdir)
        msg = ["Running...", "\t%s" % self._cmd]
        for m in msg:
            self._logger.info(m)
        retcode = call(self._cmd, shell=True)
        if retcode < 0:
            self._logger.error("join_paired_ends.py terminated by " +
                               "signal %s" % -retcode)
            sys.exit(1)
        else:
            self._logger.info("join_paired_ends.py returned %s" % retcode)
        return self._outfname

    def __build_cmd(self, infnames, outdir):
        """Build a command-line for join_paired_ends.py"""
        f1, f2 = tuple(infnames)
        cmd = ["join_paired_ends.py",
               "-f", f1,
               "-r", f2,
               "-o", outdir]
        self._cmd = ' '.join(cmd)
        self._outfname = os.path.join(outdir, "fastqjoin.join.fastq")
