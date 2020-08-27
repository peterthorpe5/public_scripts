#!/usr/bin/env python
import os


def test_if_id_is_metazoan(tax_id_of_interst,
                           final_tx_id_to_identify_up_to,
                           tax_to_filter_out):
    """function to get a list of tax id of interest from the
    tax_dictionary which is produced in the parse_function
    (parse_NCBI_nodes_tab_file)
    nodes.dmp file. . The tax id
    are added to a list for later use
    """
    # print "filtering up to =", final_tx_id_to_identify_up_to
    # print "filtering out = ", tax_to_filter_out
    if tax_id_of_interst == "N/A":
        raise ValueError("N/A as taxonomy ID")
    if tax_id_of_interst == "0":
        tax_id_of_interst == "32644"  # assign an unknown id
        return "In_filter_out_tax_id"
    if tax_id_of_interst == "7070":
        return "In_filter_out_tax_id"
    # call the function to parse nodes file and assign to variable
    # tax_dictionary = parse_NCBI_nodes_tab_file(nodes_dmp)
    # empty list to add tax id to
    # list_of_tx_id_identified_that_we_want = []
    # get the "master" parent id
    parent = tax_dictionary[tax_id_of_interst]
    # print parent
    # list_of_tx_id_identified_that_we_want.append(parent)
    while True:
    # for keys in tax_dictionary:
        # print "parent = ", parent, "\n"
        parent = tax_dictionary[parent]
        if tax_id_of_interst == "N/A":
            raise ValueError("N/A as taxonomy ID")
        # list_of_tx_id_identified_that_we_want.append(parent)
        # print list_of_tx_id_identified_that_we_want

        # 32630 is a synthetic organism
        if parent == "32630":  # 32630
            return "In_filter_out_tax_id"
            break
        if parent == "7070":
            return "In_filter_out_tax_id"
        if parent == tax_to_filter_out:
            # print "filtering out"
            return "In_filter_out_tax_id"
            break
        if parent == final_tx_id_to_identify_up_to:
            # print "......................... im here"
            return True
        elif parent == "1":
            # print "Reached the root of the tree"
            return False


