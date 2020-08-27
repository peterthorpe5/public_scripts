#!/usr/bin/env python
import os

def parse_NCBI_nodes_tab_file(folder):
    """this is a function to open nodes.dmp from the NCBI taxonomy
    database and find the parent child relationship....returns a
    dictionary for later use"""
    # open file - read.
    # nodes.dmp - this file is separated by \t|\t
    # tax_dmp_database = open(nodes_dmp, "r")
    # empty dictionary to add to parent and child (keys,vals) to
    tax_dictionary = {}
    # nodes.dmp files goes: child, parent, etc
    # merged.dmp file goes: old, new
    # In both cases, can take key as column 0 and value as column 1
    for filename in ["nodes.dmp", "merged.dmp"]:
        with open(os.path.join(folder, filename)) as handle:
            for line in handle:
                tax_info = line.replace("\n", "\t").split("\t|\t")
                # first element
                parent = tax_info[1]
                # second element
                child = tax_info[0]
                # add these to the dictionary {parent:child}
                tax_dictionary[child]= parent
    # print tax_dictionary
    # f_out = open("dictionary_test.out", "w")
    # print >> f_out, tax_dictionary
    # f_out.close()
    # print(tax_dictionary)
    return tax_dictionary


