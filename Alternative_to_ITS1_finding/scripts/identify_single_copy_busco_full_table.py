#!/usr/bin/env python
# author: Peter Thorpe and Leighton Pritchard. The James Hutton Insitute,Dundee,UK.
# 2017 Feb

# Title:
# function to identify the single copy busco genes from the full table output.
# Busco Version 1.1b

""" Example: data in:
#BUSCO_group	Status	Scaffold	Start	End	Bitscore	Length
BUSCOaEOG7ZQ0ZR	Duplicated	JOTR01000471.1	15224	43397	1902.3	1739
BUSCOaEOG71CT16	Duplicated	JOTR01000058.1	443982	461304	261.9	464
BUSCOaEOG7C8TZQ	Complete	JOTR01000184.1	422466	438678	181.3	168
BUSCOaEOG7HTVX9	Complete	JOTR01000363.1	209433	225546	192.0	318
BUSCOaEOG7BKR80	Complete	JOTR01000169.1	192730	258694	895.4	372
"""
###########################################################################


def parse_tab_outfile(busco):
    """read in the busco tab file. Reads whole file into memory.
    returns a list, one list item per busco hit.
    """
    with open(busco) as file:
        return file.read().split("\n")


def complete_single_busco(busco, genome, genom_scaff_start_stop):
    """function to split the lines and return them for
    another function to use"""
    # call the function
    duplic_frag_list = []
    busco_hits = parse_tab_outfile(busco)
    # list to append the results to
    complete_busco = []
    name_list = []

    for i in busco_hits:
        if i.startswith("#"): #allows line to have comment.
            continue
        if not i.strip():
            continue #if the last line is blank
        if len(i.split("\t")) < 7:
            # e.g. BUSCOaEOG78DM0G	Missing
            continue
        BUSCO_group, Status, Scaff, start, stop, Bitscore, Length = i.split("\t")
        # It seems that the len is always much shorter!!!!!!
        # assert Length == int(stop) - int(start), (
           #  "len = %s, star stop: %s %s, stop -start  = %d" % (Length, start, stop,
                                                           # (int(stop) - int(start))))
        if Status != "Complete":
            duplic_frag_list.append(BUSCO_group)
            continue
        if BUSCO_group  in duplic_frag_list:
            print ("you should not see this")
            continue
        if int(start) == 0:
            start = "1"
        complete_busco_line = ("\t".join([BUSCO_group, Status, Scaff,
                                        start, stop, Bitscore, Length]))
        assert BUSCO_group not in name_list
        name_list.append(BUSCO_group)
        complete_busco.append(complete_busco_line)
        # add the data to a default dictionary
        genom_scaff_start_stop[genome].append("\t".join([BUSCO_group,
                                                       Scaff,
                                                       start,
                                                       stop]))
    return complete_busco, name_list, genom_scaff_start_stop

