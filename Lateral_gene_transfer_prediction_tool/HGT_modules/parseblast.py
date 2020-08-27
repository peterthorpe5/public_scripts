#!/usr/bin/env python
import os

def parse_blast_line(blast_line_as_list, tax_coloumn):
    """function to return the bits of info which we want
    from the blast line(passed to function as a list)"""
    if len(blast_line_as_list) == 1:
        blast_line_as_list = blast_line_as_list[0]
    blast_line = blast_line_as_list
    try:
        Evalue = float(blast_line[10])
        bit_score = float(blast_line[11])
    except ValueError:
        print("\n\ncannot convert to float", blast_line[10],
               blast_line[11], "\n")
    # tax id can have a whole load of value e.g.
    # 5141;367110;510951;510952;771870.
    # Thefore we split it and take the first one
    query_name = blast_line[0]
    query_name = query_name.split("gene=")[0]
    # if "g3392" in query_name:
        # print "RAW BLAST line = ", query_name, blast_line
    percentage_identity = blast_line[2]
    description = blast_line[12]
    tax_id = blast_line[tax_coloumn].split(";")[0]
    species_sci = blast_line[-3]
    species_common = blast_line[-2]
    kingdom = blast_line[-1]
    return query_name, percentage_identity, Evalue, bit_score, \
           description, tax_id, species_sci, species_common,\
           kingdom


