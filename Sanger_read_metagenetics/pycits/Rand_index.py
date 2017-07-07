#!/usr/bin/env python
# coding: utf-8
# Title:
# script to tests rand index

# author: Peter Thorpe and Leighton Pritchard
# November 2016. The James Hutton Insitute, Dundee, UK.

import os
from sklearn.metrics.cluster import adjusted_rand_score


def open_parse(clustr):
    """function: opens file and return
    list split on \n from the file"""
    with open(clustr) as file:
        data = file.read().split("\n")
    return data


def return_rand_cluster_member(line):
    """function: split the line """
    element = line.split("\t")
    cluster = element[1]
    return cluster.rstrip()


def prepare_rand_list(infile):
    """function to repare lists of the wanted info"""
    data = open_parse(infile)
    cluster_list = []
    for line in data:
        if not line.strip():
            continue  # if the last line is blank
        if line.startswith("#"):  # dont want comment lines
            continue
        cluster = return_rand_cluster_member(line)
        cluster_list.append(cluster)
    return cluster_list


def Rand_index_cal(infile, infile2, prefix):
    """function to calcutae the rand index between
    clustering programs/ Call other functions to
    open parse the file, and return a list of results.
    requires:
    import sklearn
    from sklearn.metrics.cluster +
    import adjusted_rand_score"""
    cluster_list = prepare_rand_list(infile)
    cluster_list2 = prepare_rand_list(infile2)
    rant_result = adjusted_rand_score(cluster_list,
                                      cluster_list2)
    result = ("%s\tadjusted_rand_score =\t%f\n" %
              (prefix, rant_result))
    return result


def pairwise_comparison_Rand(file_list, outfile):
    """fucntion takes in a list of fiile names and call
    the Rand_index_cal function to compare the
    lists in the files.
    The files have format:
    read\tclusternumber
    Writes out the results to a file"""
    f_out = open(outfile, "w")
    for entry in sorted(file_list):
        for another in file_list:
            results = Rand_index_cal(entry, another,
                                     os.path.split(entry)[-1] +
                                     "_Vs_" +
                                     os.path.split(another)[-1])
            f_out.write(results)
    f_out.close()
    return results
