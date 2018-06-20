#!/usr/bin/env python

"""Test wrapper to run metapy.
The defaults are the actual test data. So it only needs to be called
as metapy.py
"""

import os
import subprocess
# using the swarm module to test for execution
from pycits import swarm
from pycits.tools import NotExecutableError
from nose.tools import nottest, assert_equal


def test_metapy():
    """metapy instantiates with cmd-line if metapy.py is in $PATH"""
    # using swarm module here just to check it is present
    cluster = swarm.Swarm("metapy.py")


def test_metapy_exec_notexist():
    """Error thrown if metapy.py executable does not exist"""
    try:
        # using swarm module here just to check it is present
        cluster = swarm.Swarm(os.path.join(".", "metapy.py"))
    except NotExecutableError:
        return True
    else:
        return False


def test_metapy_exec():
    """Run metapy.py on test data compare  to pre target.
    """
    # The default option are the actual
    # test data. So we only need to call the program
    cmd_s = "./metapy.py"
    pipe = subprocess.run(cmd_s, shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          check=True)
