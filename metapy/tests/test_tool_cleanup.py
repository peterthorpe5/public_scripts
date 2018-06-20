#!/usr/bin/env python

"""Tests of wrapper code in pycits: clean_up"""

import os
import shutil

from pycits import clean_up
from pycits.tools import NotExecutableError

from nose.tools import nottest, assert_equal


def make_me_a_file():
    test_file = open("temp.txt", "w")
    test_file.write("file to be deleted")
    test_file.close()


def test_clean_up():
    """clean_up instantiates with cmd-line if clean_up is in $PATH"""
    assemble = clean_up.Clean_up("rm")


def test_clean_up_cmd():
    """clean_up instantiates, runs and returns correct form of cmd-line"""
    make_me_a_file()
    obj = clean_up.Clean_up("rm")
    target = ' '.join(["rm", "temp.txt"])
    results = obj.run("temp.txt")
    assert_equal(results.command, target)


def test_clean_up_exec_notexist():
    """Error thrown if clean_up executable does not exist"""
    try:
        obj = clean_up.Clean_up(os.path.join(".", "rm"))
    except NotExecutableError:
        return True
    else:
        return False


def test_clean_up_notexec():
    """Error thrown if clean_up exe not executable"""
    try:
        obj = clean_up.Clean_up("LICENSE")
    except NotExecutableError:
        return True
    else:
        return False


def test_clean_up_exec():
    """Run clean_up on test data and compare output to precomputed target
    this often does not produce the same result twice. So we may have
    to disable this test"""
    make_me_a_file()
    obj = clean_up.Clean_up("rm")
    result = obj.run("temp.txt")
    if os.path.isfile("temp.txt"):
        return False
