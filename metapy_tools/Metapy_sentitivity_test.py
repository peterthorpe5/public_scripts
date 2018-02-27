# script to quantify sensitivity and other metrics for metapy output
# paste the results from the output of interest into the variable called:
# FOUND_STR
# Also, make sure you test the FOUND_STR against the correct TEST_LIST!!

DNAMIXUNDIL = """capsici
obscura
castaneae
siskiyouensis
plurivora
foliorum
rubi
fallax
cactorum
boehmeriae""".split()

DNAMIX2 = """pseudosyringae
ilicis
austrocedri 
kernoviae
cactorum
chlamydospora
gonapodyides
obscura
ramorum
cinnamomi
lateralis
cambivora
syringae
plurivora
boehmeriae""".split()

DNAMIX1 = """idaei
capsici
plurivora
palmivora
castaneae
megasperma
rubi
cryptogea
fallax
boehmeriae""".split()

PCRMIX = """idaei
capsici
plurivora
palmivora
castaneae
megasperma
rubi
cryptogea
fallax
boehmeriae""".split()


#######################################################
# CHANGE YOUR COMPARISON TEST LIST HERE
# DNAMIXUNDIL, DNAMIX2, DNAMIX1, PCRMIX
TEST_LIST = DNAMIX1

# Paste you results in here
FOUND_STR = """1	7Phytophthora_uliginosa_P10413 7Phytophthora_cinnamomi_var_robiniae_P16351_P.asiatica_P16351_P.europaea_BR1072 7Phytophthora_melonis_CBS58269_7Phytophthora_sinensis_CBS55788 7Phytophthora_niederhauserii_P10279 7Phytophthora_sojae_CBS38261 Phytophthora_cinnamomi_TW102 7Phytophthora_cinnamomi_CBS14422 6Phytophthora_chlamydospora_VHS6595_EU301159 6Phytophthora_sp_nov_CBS114338_Phytophthora_taxon_PgChlamydo_VHS3753_EU301160 6Phytophthora_gonapodyides_voucher_HQ643236_1 6Phytophthora_mississippiae_ATCC_MYA_4946 Phytophthora_gibbosa_CBS127951_HQ012933_Phytophthora_gregata_CBS127952_HQ012942 6Phytophthora_megasperma_IMI133317_AF266794_P.crassamura_DDS3432 6Phytophthora_sp_nov_BR333 6Phytophthora_taxon_hungarica_JA566 Phytophthora_taxon_raspberry_P1050_AF541905 P._borealis_AKWA58.10708_HM004232 Phytophthora_x_stagnum_clone_type1 Phytophthora_x_stagnum_clone_type2 1Phytophthora_phaseoli_CBS55688 1Phytophthora_versiformis_MJ5 5Phytophthora_katsurae_CBS58785 5Phytophthora_cocois_ICMP16948__KP295304 5Phytophthora_sp_novaeguinea_P1256_5P.agathidicida_KP295308_5P.castaneae_TW270 4Phytophthora_alticola_CMW34279_HQ013214_4P.arenaria__CBS125800_HQ13215 4Phytophthora_palmivora_CBS17926_4Phytophthora_arecae_CBS30562 8Phytophthora_himalayensis_CBS35759_8Phytophthora_erythroseptica_CBS111343 8Phytophthora_primulae_CBS110162_8Phytophthora_aff_primulae_P6817 2Phytophthora_capsici_CBS12823_2Phytophthora_sp_glovera_P10618 2Phytophthora_mexicana_P0646	1224
2	8Phytophthora_sp_kelmania_CBS30762	364
3	2Phytophthora_citricola_CBS22188 2Phytophthora_pini_CBS18125 2Phytophthora_plurivora_CBS37961 2Phytophthora_acerina_B035 2Phytophthora_pachypleura_RHS2474 2Phytophthora_sp._FFM22.2	301
4	1Phytophthora_pseudotsugae01_voucherCBS44484 1Phytophthora_cactorum_CBS113344 1Phytophthora_idaei_P6767 1Phytophthora_pseudotsugae02_voucherP10339	301
5	7Phytophthora_fragariae_CBS20946 7Phytophthora_x_multiformis_BBA_PAM54KJ755097 7Phytophthora_rubi_CBS109892	293
6	10Phytophthora_boehmeriae_CBS29129	276
7	9Phytophthora_fallax_P10722	151
8	1Phytophthora_tabaci_CBS30529	2
10	3Phytophthora_nemorosa_CBS114870 3Phytophthora_pseudosyringae_391716	1
17	7Phytophthora_intricata_CBS141211_KU517155	1

"""

FOUND_STR = FOUND_STR.rstrip("\n")

# THIS WAS MISSPELLED IN AN OLD GENBANK ENTRY
FOUND_STR = FOUND_STR.replace("austrocedrae", "austrocedri")

FOUND_LIST = FOUND_STR.split("\n")
    
########################################################
TP = set([])
TP_list = []
FP_set_temp = set([])
seen_wanted = set([])
FN_set = set([])
FN = []
total_found_count = 0
TP_line = ""
TEST_STR = ""

for SPIKE_SPECIES in TEST_LIST:
    TEST_STR = TEST_STR + "SPIKE_SPECIES"
    # iterate through the spike test mix
    # add the names to a set to get the length later
    seen_wanted.add(SPIKE_SPECIES.strip())
    if SPIKE_SPECIES.strip() in FOUND_STR:
        TP_list.append(SPIKE_SPECIES)            
    if SPIKE_SPECIES.strip() not in FOUND_STR:
        FN.append(SPIKE_SPECIES.strip())

# setting up a str of the "found" enteries
for result in FOUND_LIST:
    for TP_found in TP_list:
        if TP_found in result:
            TP_line = TP_line + result

# convert the TP list to a set
for true_entry in TP_list:
    TP.add(true_entry)

# convert the FN list to a set    
for false_entry in FN:
    if false_entry in TP:
        continue
    FN_set.add(false_entry)

# count the total found in the experiment
for TEST_FOUND in FOUND_LIST:
    total_found_count += 1

# get a temp FN set
for found in FOUND_LIST:
    species = found.split("_")[2]
    if species not in TEST_STR:
        # make sure these are not clustered with the true positives. 
        if species not in TP_line:
            FP_set_temp.add(species)

# remove any species from the temp FP set with have TP and FN sets
FP_set = FP_set_temp.difference(TP, FN_set)
# iterate through the found list - this can be formatted differently
for found in FOUND_LIST:
    total_found_count += 1

###############################################################################
print "FALSE negative", FN_set, len(FN_set)
print "True positives = ", len(TP), TP
FP = len(FP_set)
print "FALSE postivies = ", FP_set
sens = float(len(TP))/(len(TP) + len(FN_set))
prec = float(len(TP))/(len(TP) + FP)
FDR = float(FP)/(len(TP) + FP)
FNR = float(len(FN_set))/(len(TP) + len(FN_set))
print "sensitivity = ", sens
print "Precision = ", prec
print "False discovery rate = ", FDR
print "False negative rate = ", FNR

print "THIS data set found: ", FOUND_STR.count("Phy"), "Phytophthora species"

print "%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%s Phytophthora species" %(len(TP), FP, sens,
                                             prec, FDR,
                                             FNR, FOUND_STR.count("Phy"))



                        
                        

