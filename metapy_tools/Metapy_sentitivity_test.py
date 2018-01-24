# script to quantify sensitivity and other metrics for metapy output

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
TEST_LIST = DNAMIX2

# Paste you results in here
FOUND_STR = """2	2_Phytophthora_mexicana_HQ261620_1 2_Phytophthora_capsici_AF266787_2b_P._glovera_AF279124_1	471
3	8d_Phytophthora_obscura_HQ917910_1	429
5	2_Phytophthora_siskiyouensis_EF523386_1	173
6	8c_Phytophthora_foliorum_EF120469_1	150
7	2_Phytophthora_citricola_E_AF266788_2_P._acerina_JX951285_1 2_Phytophthora_citricola_sensu_stricto_FJ237526_1 2_Phytophthora_plurivora_FJ665225_1 2_Phytophthora_pini_FJ392322_1 2_Phytophthora_pachypleura_KC855330_1 2c_Phytophthora_capensis_GU191231_1	134

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

print "%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f" %(len(TP), FP, sens, prec, FDR, FNR)



                        
                        

