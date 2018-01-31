#!/usr/bin/env python
#
# vsearch.py
#
# Wrapper code for the VSEARCH package, available from:
# https://github.com/torognes/vsearch
#
# (c) The James Hutton Institute 2016-2017
# Author: Leighton Pritchard and Peter Thorpe

"""vsearch.py

Provides the Vsearch class as a common wrapper to several Vsearch commands
used in the pycits/metapy workflow:

vsearch --cluster_fast FILENAME --id 0.97 --centroids FILENAME
vsearch --derep_fulllength FILENAME --output FILENAME
vsearch --usearch_global FILENAME --db FILENAME --id 0.97 --alnout FILENAME

At the terminal, commands are accessed by the first argument to vsearch, and
this argument controls which parameters are expected by the program.

In this implementation, we abstract out the running of any VSEARCH command
to the structure:

vsearch <MODE> <INPUTFILE> <OUTPUTFILE> <PARAMETERS>

The <MODE> is used for one of the VSEARCH command options (e.g.
--cluster_fast) and is expected as a string including the two hyphens.

With this abstraction, options are passed to an instantiated object when the
.run() method is called. <INPUTFILE> is passed in place of FILENAME (in
the VSEARCH usage examples), and <OUTPUTFILE> is used for one of the
output options - exactly which depends on the chosen <MODE>.

All other parameter are passed when calling .run() as a dictionary, keyed
by argument flag (e.g. --db) with value of the argument value (such as a
filename). All parameters for a run - whether optional or required for
VSEARCH command completion - must be passed in <PARAMETERS>.

Example usage:

cluster_params = {'--blast6out': "my_cluster.blast6",
                  '--id': 0.96,
                  '--db': "my_sequences.fasta",
                  '--threads': 2}
mode = '--usearch_global'
vsearch_exe = vsearch.Vsearch("vsearch")
result = vsearch_exe.run(mode, "my_input_sequences.fasta",
                         "my_output_cluster.uc",
                         cluster_params)
"""

import os
import subprocess
from collections import namedtuple
from .tools import is_exe, NotExecutableError

# factory class for Vsearch class returned values
Results_derep = namedtuple("Results",
                           "command outfilename stdout stderr")

Results_cluster = namedtuple("Results",
                             "command outfile_uc outfile_b6 stdout stderr")

Results_cluster_fast = namedtuple("Results",
                                  "command outfile_uc outfile_b6 outfile_msa " +
                                  "outfile_consensus centroids stdout stderr")


class VsearchError(Exception):
    """Exception raised when Vsearch fails"""
    def __init__(self, message):
        self.message = message


class Vsearch(object):
    """Class for working with VSEARCH"""
    def __init__(self, exe_path):
        """Instantiate with location of executable"""
        if not is_exe(exe_path):
            msg = "{0} is not an executable".format(exe_path)
            raise NotExecutableError(msg)
        self._exe_path = exe_path

        # Send to different command-builder depending on operation, and
        # support different output depending on operation
        self._builders = {'--derep_fulllength': self.__build_cmd_derep,
                          '--usearch_global': self.__build_cmd_cluster,
                          '--cluster_fast': self.__build_cmd_cluster_fast}
        self._returnval = {'--derep_fulllength': self.__return_derep,
                           '--usearch_global': self.__return_cluster,
                           '--cluster_fast': self.__return_cluster_fast}
        
    def run(self, mode, infile, output, params=None, dry_run=False):
        """Run VSEARCH in the prescribed mode

        - mode: one of the VSEARCH arguments for mode of operation
        - infile: path to input file
        - output: path to output directory or filename (depends on operation)
        - params: dictionary of parameters for the VSEARCH operation,
                  whether these are required for an operation may vary
        - dry_run: if True returns cmd-line but does not run

        Returns namedtuple with a form reflecting the operation requested:
          --derep_fulllength
            "command outfile stdout stderr"
        """
        self._params = params
        try:
            self._builders[mode](infile, output)
        except KeyError:
            msg = "VSEARCH mode {0} not supported".format(mode)
            raise VsearchError(msg)
        if dry_run:
            pipe = None
        else:
            pipe = subprocess.run(self._cmd, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True)
        return self._returnval[mode](pipe)

    def __build_cmd_derep(self, infile, outfname):
        """Run VSEARCH to dereplicate input sequences

        The --sizeout option is enforced, to add sequence abundance and sort
        """
        self._outfile = outfname
        cmd = ['vsearch', '--derep_fulllength', infile,
               '--output', self._outfile,
               '--sizeout']
        if self._params is not None:
            for k, v in self._params.items():
                cmd.append('--{0} {1}'.format(k, v))
        self._cmd = ' '.join(cmd)

    def __return_derep(self, pipe):
        """Returns output values for dereplication of input sequences.

        The return value is a Results_derep namedtuple
        """
        if pipe is None:  # It was a dry run
            results = Results_derep(self._cmd, self._outfile, None, None)
        else:
            results = Results_derep(self._cmd, self._outfile,
                                    pipe.stdout, pipe.stderr)
        return results

    def __build_cmd_cluster(self, infile, outfname):
        """Run VSEARCH to cluster input sequences

        The command expects the output .uc file to be specified, but all other
        options are expected in params
        """
        # Some parameters are required for operation
        for required in ['--id', '--db']:
            if required not in self._params:
                msg = "Required parameter {0} not passed".format(required)
                raise VsearchError(msg)
        self._outfile = outfname
        cmd = ['vsearch', '--usearch_global', infile,
               '--uc', self._outfile]
        if self._params is not None:
            for k, v in sorted(self._params.items()):
                cmd.append('{0} {1}'.format(k, v))
        self._cmd = ' '.join(cmd)

    def __return_cluster(self, pipe):
        """Returns output values for clustering input sequences.

        The return value is a Results_cluster namedtuple
        """
        # The blast6 output may not be defined in parameters, but is needed
        # for the return value
        if '--blast6out' not in self._params:
            self._params['--blast6out'] = None
        if pipe is None:  # It was a dry run
            results = Results_cluster(self._cmd, self._outfile,
                                      self._params['--blast6out'],
                                      None, None)
        else:
            results = Results_cluster(self._cmd, self._outfile,
                                      self._params['--blast6out'],
                                      pipe.stdout, pipe.stderr)
        return results
    

    def __build_cmd_cluster_fast(self, infile, outfname):
        """Run VSEARCH to cluster_fast input sequences

        The command expects the output .uc file to be specified, but all other
        options are expected in params
        """
        # Some parameters are required for operation
        for required in ['--id', '--centroids']:
            if required not in self._params:
                msg = "Required parameter {0} not passed".format(required)
                raise VsearchError(msg)
        self._outfile = outfname
        cmd = ['vsearch', '--cluster_fast', infile,
               '--uc', self._outfile]
        if self._params is not None:
            for k, v in sorted(self._params.items()):
                cmd.append('{0} {1}'.format(k, v))
        self._cmd = ' '.join(cmd)

    def __return_cluster_fast(self, pipe):
        """Returns output values for cluster_fast of input sequences.

        The return value is a Results_cluster namedtuple
        """
        # Optional output may not be defined in parameters, but is needed
        # for the return value
        for param in ['--blast6out', '--msaout', '--consout']:
            if param not in self._params:
                self._params[param] = None
        if pipe is None:  # It was a dry run
            results = Results_cluster_fast(self._cmd, self._outfile,
                                           self._params['--blast6out'],
                                           self._params['--msaout'],
                                           self._params['--consout'],
                                           self._params['--centroids'],
                                           None, None)
        else:
            results = Results_cluster_fast(self._cmd, self._outfile,
                                           self._params['--blast6out'],
                                           self._params['--msaout'],
                                           self._params['--consout'],
                                           self._params['--centroids'],
                                           pipe.stdout, pipe.stderr)
        return results
