#!/usr/bin/env python3
# title: compare what species were found using sets
# in the cluster
# author: Peter Thorpe and Leighton Pritchard
# September 2016. The James Hutton Insitute, Dundee, UK.

# imports
import os
import sys
from optparse import OptionParser
import os
from sys import stdin,argv
from collections import defaultdict
###############################################################


def stop_err( msg ):
    """stop function"""
    sys.stderr.write("%s\n" % msg)
    sys.exit(1)


def open_parse(infile):
    """function: opens file and return
    list split on \n from the file"""
    with open(infile) as file:
        data = file.read().split("\n")
    return data


def parse_line(line):
    """function to parse line"""
    if not line.strip():
        return 0 #if the last line is blank
    if line.startswith("#"):
        return 0
    if line.startswith("    "):
        return 0
    line = line.rstrip()
    out_put_str = ""
    result_line = line.rstrip("\n").split("\t")
    cluster_num = result_line[0]
    species = result_line[1]
    reads = result_line[2]
    return cluster_num, species, reads


def parse_list(in_list, out_file):
    """script to open up a tab separeted clustering output
    and identify the species in the clustering"""
    spe_fnd_dic = defaultdict(set)
    in_list = open_parse(in_list)
    for entry in in_list:
        if entry == "":
            continue
        program, result_file = entry.split("\t")
        data = open_parse(result_file)
        if data == 0:
            continue
        for line in data:
            if parse_line(line) == 0:
                continue
            else:
                cluster_num, species, reads = parse_line(line)
                species = species.split()
                for single in species:
                    if "_abundance=1" in single:
                        single = single.split("_abundance=1")[0]
                    spe_fnd_dic[program].add(single)
    return spe_fnd_dic


def compare_found(spe_fnd_dic, out_file):
    """function to compare a disctionary of sets to identify what was,
    found and what was unique to each tool"""
    number_of_items = 0
    names = []
    f_out = open(out_file, "w")
    common_in_all = spe_fnd_dic["bowtie"] & \
                    spe_fnd_dic["swarm"] \
                    & spe_fnd_dic["cdhit"] \
                    & spe_fnd_dic["vsearch"] & \
                    spe_fnd_dic["vse_faclus"] \
                    & spe_fnd_dic["blastclust"]
    out = "#%d Species found in all methods:\t%s\n\n" % (len(common_in_all),
                                                         str(common_in_all))
    f_out.write(out)
    for program, species in spe_fnd_dic.items():
        number_of_items += 1
        names.append(program)
        out_str = "%s method found %d species:\t%s" % (program, len(species),
                                                       species)
        f_out.write(out_str)
        f_out.write("\n")
    f_out.write("\n#comparison of 4 methods at a time\n")
    for n1 in names:
        for n2 in names:
            if n1 == n2:
                continue
            for n3 in names:
                if n3 == n2:
                    continue
                if n3 == n1:
                    continue
                for n4 in names:
                    if n4 == n1:
                        continue
                    if n4 == n2:
                        continue
                    if n4 == n3:
                        continue
                    common_var = spe_fnd_dic[n1].intersection(spe_fnd_dic[n2],
                                                              spe_fnd_dic[n3],
                                                              spe_fnd_dic[n4])
                    diff_var = spe_fnd_dic[n1].difference(spe_fnd_dic[n2],
                                                          spe_fnd_dic[n3],
                                                          spe_fnd_dic[n4])
                    fmd = "%s_vs_%s_%s_%s_%d_COMMON:\t%s\n" %(n1, n2,
                                                              n3, n4,
                                                              len(common_var),
                                                              str(common_var))
                    f_out.write(fmd)
                    fmd = "%s_vs_%s_%s_%s_%d_DIFFERENCE:\t%s\n" %(n1, n2,
                                                                  n3, n4,
                                                                  len(diff_var),
                                                                  str(diff_var))
                    f_out.write(fmd)


    #####################################
    f_out.write("\n#comparison of 3 methods at a time\n")
    for n1 in names:
        for n2 in names:
            if n1 == n2:
                continue
            for n3 in names:
                if n3 == n2:
                    continue
                if n3 == n1:
                    continue
                common_var = spe_fnd_dic[n1].intersection(spe_fnd_dic[n2],
                                                          spe_fnd_dic[n3])
                diff_var = spe_fnd_dic[n1].difference(spe_fnd_dic[n2],
                                                      spe_fnd_dic[n3])
                formmatted = "%s_vs_%s_%s_%d_COMMON:\t%s\n" %(n1, n2, n3,
                                                              len(common_var),
                                                              str(common_var))
                f_out.write(formmatted)
                formated = "%s_vs_%s_%s_%d_DIFFERENCE:\t%s\n" %(n1, n2, n3,
                                                                len(diff_var),
                                                                str(diff_var))
                f_out.write(formated)
    ###########################
    f_out.write("\n#comparison of 2 methods at a time\n")
    for n1 in names:
        for n2 in names:
            if n1 == n2:
                continue
            common_var = spe_fnd_dic[n1].intersection(spe_fnd_dic[n2])
            diff_var = spe_fnd_dic[n1].difference(spe_fnd_dic[n2])
            formmatted = "%s_vs_%s_%d_COMMON:\t%s\n" %(n1, n2,
                                                       len(common_var),
                                                       str(common_var))
            f_out.write(formmatted)
            formated = "%s_vs_%s_%d_DIFFERENCE:\t%s\n" %(n1, n2,
                                                         len(diff_var),
                                                         str(diff_var))
            f_out.write(formated)
    f_out.close()

#############################################################################
# to run the script

usage = """usage :

script takes strings. e.g.

'swarm\t/PATHTO/DNAMIX_S95_L001_Swarm_d1/swarm_results.RESULTS\n
cdhit\t/PATHTO/DNAMIX_S95_L001_cd_hit_0.99/cdhit_results_0.99.RESULTS\n'

name of tools \t path to restuls file.

The script will compare what species were found in each tools and
report common, unique speices found per tools

example of result file format:

# cluster_number	species	number_of_reads_hitting_species
1	2Phytophthora_capsici_CBS12823	22
2	2Phytophthora_plurivora_CBS37961	141
3	2Phytophthora_sp_glovera_P10618	14
4	4Phytophthora_arecae_CBS30562	1
5	5Phytophthora_sp_novaeguinea_P1256	46
6	7Phytophthora_rubi_CBS109892	102
7	Phytophthora_megasperma_IMI133317_AF266794	22
    # Fasta file assembled seq summary:
    # min_contig = 113 max_contig = 233 avg_contig = 189
    # Total number of assembled sequences = 3414
    # number of assembled-reads clustering with database = 348
    # number of starting reads = 5159
    # percent of assembled-reads clustering with database = 10.19
    # db_fa_summary:
    # dbmin_contig = 159 dbmax_contig = 224 dbavg_contig = 195
    # number of db_seq = 153

"""

parser = OptionParser(usage=usage)

parser.add_option("--in_list", dest="in_list",
                  default=None,
                  help="list strings of files: ['program\tfile']")

parser.add_option("-o", "--out_prefix", dest="out_file",
                  default="summarise_clusters.out",
                  help="prefix to the output filenames")

(options, args) = parser.parse_args()
# --in_list
in_list = options.in_list
# -o
out_file = options.out_file


###################################################################
# run the program
spe_fnd_dic = parse_list(in_list, out_file)
compare_found(spe_fnd_dic, out_file)
=======
#!/usr/bin/env python
# title: compare what species were found using sets
# in the cluster
# author: Peter Thorpe and Leighton Pritchard
# September 2016. The James Hutton Insitute, Dundee, UK.

# imports
import os
import sys
from optparse import OptionParser
import os
from sys import stdin,argv
from collections import defaultdict
###############################################################


def stop_err( msg ):
    """stop function"""
    sys.stderr.write("%s\n" % msg)
    sys.exit(1)


def open_parse(infile):
    """function: opens file and return
    list split on \n from the file"""
    with open(infile) as file:
        data = file.read().split("\n")
    return data


def parse_line(line):
    """function to parse line"""
    if not line.strip():
        return 0 #if the last line is blank
    if line.startswith("#"):
        return 0
    if line.startswith("    "):
        return 0
    line = line.rstrip()
    out_put_str = ""
    result_line = line.rstrip("\n").split("\t")
    cluster_num = result_line[0]
    species = result_line[1]
    reads = result_line[2]
    return cluster_num, species, reads


def parse_list(in_list, out_file):
    """script to open up a tab separeted clustering output
    and identify the species in the clustering"""
    spe_fnd_dic = defaultdict(set)
    in_list = open_parse(in_list)
    for entry in in_list:
        if entry == "":
            continue
        program, result_file = entry.split("\t")
        data = open_parse(result_file)
        if data == 0:
            continue
        for line in data:
            if parse_line(line) == 0:
                continue
            else:
                cluster_num, species, reads = parse_line(line)
                species = species.split()
                for single in species:
                    if "_abundance=1" in single:
                        single = single.split("_abundance=1")[0]
                    spe_fnd_dic[program].add(single)
    return spe_fnd_dic


def compare_found(spe_fnd_dic, out_file):
    """function to compare a disctionary of sets to identify what was,
    found and what was unique to each tool"""
    number_of_items = 0
    names = []
    f_out = open(out_file, "w")
    common_in_all = spe_fnd_dic["bowtie"] & \
                    spe_fnd_dic["swarm"] \
                    & spe_fnd_dic["cdhit"] \
                    & spe_fnd_dic["vsearch"] & \
                    spe_fnd_dic["vse_faclus"] \
                    & spe_fnd_dic["blastclust"]
    out = "#%d Species found in all methods:\t%s\n\n" % (len(common_in_all),
                                                         str(common_in_all))
    f_out.write(out)
    for program, species in spe_fnd_dic.items():
        number_of_items += 1
        names.append(program)
        out_str = "%s method found %d species:\t%s" % (program, len(species),
                                                       species)
        f_out.write(out_str)
        f_out.write("\n")
    f_out.write("\n#comparison of 4 methods at a time\n")
    for n1 in names:
        for n2 in names:
            if n1 == n2:
                continue
            for n3 in names:
                if n3 == n2:
                    continue
                if n3 == n1:
                    continue
                for n4 in names:
                    if n4 == n1:
                        continue
                    if n4 == n2:
                        continue
                    if n4 == n3:
                        continue
                    common_var = spe_fnd_dic[n1].intersection(spe_fnd_dic[n2],
                                                              spe_fnd_dic[n3],
                                                              spe_fnd_dic[n4])
                    diff_var = spe_fnd_dic[n1].difference(spe_fnd_dic[n2],
                                                          spe_fnd_dic[n3],
                                                          spe_fnd_dic[n4])
                    fmd = "%s_vs_%s_%s_%s_%d_COMMON:\t%s\n" %(n1, n2,
                                                              n3, n4,
                                                              len(common_var),
                                                              str(common_var))
                    f_out.write(fmd)
                    fmd = "%s_vs_%s_%s_%s_%d_DIFFERENCE:\t%s\n" %(n1, n2,
                                                                  n3, n4,
                                                                  len(diff_var),
                                                                  str(diff_var))
                    f_out.write(fmd)


    #####################################
    f_out.write("\n#comparison of 3 methods at a time\n")
    for n1 in names:
        for n2 in names:
            if n1 == n2:
                continue
            for n3 in names:
                if n3 == n2:
                    continue
                if n3 == n1:
                    continue
                common_var = spe_fnd_dic[n1].intersection(spe_fnd_dic[n2],
                                                          spe_fnd_dic[n3])
                diff_var = spe_fnd_dic[n1].difference(spe_fnd_dic[n2],
                                                      spe_fnd_dic[n3])
                formmatted = "%s_vs_%s_%s_%d_COMMON:\t%s\n" %(n1, n2, n3,
                                                              len(common_var),
                                                              str(common_var))
                f_out.write(formmatted)
                formated = "%s_vs_%s_%s_%d_DIFFERENCE:\t%s\n" %(n1, n2, n3,
                                                                len(diff_var),
                                                                str(diff_var))
                f_out.write(formated)
    ###########################
    f_out.write("\n#comparison of 2 methods at a time\n")
    for n1 in names:
        for n2 in names:
            if n1 == n2:
                continue
            common_var = spe_fnd_dic[n1].intersection(spe_fnd_dic[n2])
            diff_var = spe_fnd_dic[n1].difference(spe_fnd_dic[n2])
            formmatted = "%s_vs_%s_%d_COMMON:\t%s\n" %(n1, n2,
                                                       len(common_var),
                                                       str(common_var))
            f_out.write(formmatted)
            formated = "%s_vs_%s_%d_DIFFERENCE:\t%s\n" %(n1, n2,
                                                         len(diff_var),
                                                         str(diff_var))
            f_out.write(formated)
    f_out.close()

#############################################################################
# to run the script

usage = """usage :

script takes strings. e.g.

'swarm\t/PATHTO/DNAMIX_S95_L001_Swarm_d1/swarm_results.RESULTS\n
cdhit\t/PATHTO/DNAMIX_S95_L001_cd_hit_0.99/cdhit_results_0.99.RESULTS\n'

name of tools \t path to restuls file.

The script will compare what species were found in each tools and
report common, unique speices found per tools

example of result file format:

# cluster_number	species	number_of_reads_hitting_species
1	2Phytophthora_capsici_CBS12823	22
2	2Phytophthora_plurivora_CBS37961	141
3	2Phytophthora_sp_glovera_P10618	14
4	4Phytophthora_arecae_CBS30562	1
5	5Phytophthora_sp_novaeguinea_P1256	46
6	7Phytophthora_rubi_CBS109892	102
7	Phytophthora_megasperma_IMI133317_AF266794	22
    # Fasta file assembled seq summary:
    # min_contig = 113 max_contig = 233 avg_contig = 189
    # Total number of assembled sequences = 3414
    # number of assembled-reads clustering with database = 348
    # number of starting reads = 5159
    # percent of assembled-reads clustering with database = 10.19
    # db_fa_summary:
    # dbmin_contig = 159 dbmax_contig = 224 dbavg_contig = 195
    # number of db_seq = 153

"""

parser = OptionParser(usage=usage)

parser.add_option("--in_list", dest="in_list",
                  default=None,
                  help="list strings of files: ['program\tfile']")

parser.add_option("-o", "--out_prefix", dest="out_file",
                  default="summarise_clusters.out",
                  help="prefix to the output filenames")

(options, args) = parser.parse_args()
# --in_list
in_list = options.in_list
# -o
out_file = options.out_file


###################################################################
# run the program
spe_fnd_dic = parse_list(in_list, out_file)
compare_found(spe_fnd_dic, out_file)
>>>>>>> b5a3cd9f461b5e65b0b0748712c44268bd2236a1
