cd $HOME/scratch/tree_health/Illumina_results_2017Nov_first_run/

#cp /mnt/shared/projects/Phytophthora/201711_ITS1_metabarcoding_1_THAPBI/static/RawData/* ./

source ~/misc_python/THAPBI/THAPBI-pycits/venv-THAPBI-pycits-py3/bin/activate
python generate_metapy.py
##########################################################################################
wait
python interogate_results.py Control_C1
python interogate_results.py Control_C2
python interogate_results.py Control_C3
python interogate_results.py Control_C4
python interogate_results.py Control

python interogate_results.py austrocedri
python interogate_results.py austrocedrae
cat austrocedriidentified.out austrocedraeidentified.out > austrocedrae_identified.out
rm austrocedriidentified.out austrocedraeidentified.out
python interogate_results.py ramorum
python interogate_results.py lateralis
python interogate_results.py cactorum
python interogate_results.py austrocedri
python interogate_results.py gonapodyides

##########################################################################################

#mkdir reads

mv *fastq* ./reads

#  dont delete for now
#rm -rf reads

mkdir swarm_results

mkdir bowtie_results

############################################################################################

cp ./1_Phytophthora_cactorum_AF266772identified.out/1_Phytophthora_cactorum_AF266772identified.out_Swarm_d1/1_Phytophthora_cactorum_AF266772identified.out_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GL6-0x_S72_L001_RESULTS/GL6-0x_S72_L001_Swarm_d1/GL6-0x_S72_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N012-161201-5A-B_S51_L001_RESULTS/N012-161201-5A-B_S51_L001_Swarm_d1/N012-161201-5A-B_S51_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N06-161020-7B-F_S19_L001_RESULTS/N06-161020-7B-F_S19_L001_Swarm_d1/N06-161020-7B-F_S19_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./8c_Phytophthora_lateralis_AF266804identified.out/8c_Phytophthora_lateralis_AF266804identified.out_Swarm_d1/8c_Phytophthora_lateralis_AF266804identified.out_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GL6-100x_S96_L001_RESULTS/GL6-100x_S96_L001_Swarm_d1/GL6-100x_S96_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N012-161201-5B-B_S63_L001_RESULTS/N012-161201-5B-B_S63_L001_Swarm_d1/N012-161201-5B-B_S63_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N07-161028-1B-F_S31_L001_RESULTS/N07-161028-1B-F_S31_L001_Swarm_d1/N07-161028-1B-F_S31_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./8c_Phytophthora_ramorum_JQ653034identified.out/8c_Phytophthora_ramorum_JQ653034identified.out_Swarm_d1/8c_Phytophthora_ramorum_JQ653034identified.out_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GL6-10x_S84_L001_RESULTS/GL6-10x_S84_L001_Swarm_d1/GL6-10x_S84_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N012-161201-6B-B_S75_L001_RESULTS/N012-161201-6B-B_S75_L001_Swarm_d1/N012-161201-6B-B_S75_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N07-161028-3B-F_S55_L001_RESULTS/N07-161028-3B-F_S55_L001_Swarm_d1/N07-161028-3B-F_S55_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./8d_Phytophthora_austrocedri_DQ995184identified.out/8d_Phytophthora_austrocedri_DQ995184identified.out_Swarm_d1/8d_Phytophthora_austrocedri_DQ995184identified.out_swarm_results_1.RESULTS ./swarm_results/
 
cp ./interogate_results.py/interogate_results.py_Swarm_d1/interogate_results.py_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N012-161201-6C-B_S87_L001_RESULTS/N012-161201-6C-B_S87_L001_Swarm_d1/N012-161201-6C-B_S87_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N07-161028-3C-B_S67_L001_RESULTS/N07-161028-3C-B_S67_L001_Swarm_d1/N07-161028-3C-B_S67_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./all_results_files.txt/all_results_files.txt_Swarm_d1/all_results_files.txt_swarm_results_1.RESULTS ./swarm_results/
 
cp ./log_files/log_files_Swarm_d1/log_files_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N013-161205-10B-B_S76_L001_RESULTS/N013-161205-10B-B_S76_L001_Swarm_d1/N013-161205-10B-B_S76_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N07-161028-3C-F_S79_L001_RESULTS/N07-161028-3C-F_S79_L001_Swarm_d1/N07-161028-3C-F_S79_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./bowtie_results/bowtie_results_Swarm_d1/bowtie_results_swarm_results_1.RESULTS ./swarm_results/
 
cp ./metapy.sh/metapy.sh_Swarm_d1/metapy.sh_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N013-161205-1C-B_S4_L001_RESULTS/N013-161205-1C-B_S4_L001_Swarm_d1/N013-161205-1C-B_S4_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N07-161028-5A-F_S91_L001_RESULTS/N07-161028-5A-F_S91_L001_Swarm_d1/N07-161028-5A-F_S91_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./Control_C1identified.out/Control_C1identified.out_Swarm_d1/Control_C1identified.out_swarm_results_1.RESULTS ./swarm_results/
 
cp ./mv.sh/mv.sh_Swarm_d1/mv.sh_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N013-161205-3A-B_S16_L001_RESULTS/N013-161205-3A-B_S16_L001_Swarm_d1/N013-161205-3A-B_S16_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N07-161028-5C-F_S8_L001_RESULTS/N07-161028-5C-F_S8_L001_Swarm_d1/N07-161028-5C-F_S8_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./Control_C2identified.out/Control_C2identified.out_Swarm_d1/Control_C2identified.out_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N010-161128-2A-F_S38_L001_RESULTS/N010-161128-2A-F_S38_L001_Swarm_d1/N010-161128-2A-F_S38_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N013-161205-5A-B_S28_L001_RESULTS/N013-161205-5A-B_S28_L001_Swarm_d1/N013-161205-5A-B_S28_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-1A-B_S32_L001_RESULTS/N09-161117-1A-B_S32_L001_Swarm_d1/N09-161117-1A-B_S32_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./Control_C3identified.out/Control_C3identified.out_Swarm_d1/Control_C3identified.out_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N010-161128-4A-F_S50_L001_RESULTS/N010-161128-4A-F_S50_L001_Swarm_d1/N010-161128-4A-F_S50_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N013-161205-5B-B_S40_L001_RESULTS/N013-161205-5B-B_S40_L001_Swarm_d1/N013-161205-5B-B_S40_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-1A-F_S44_L001_RESULTS/N09-161117-1A-F_S44_L001_Swarm_d1/N09-161117-1A-F_S44_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./Control_C4identified.out/Control_C4identified.out_Swarm_d1/Control_C4identified.out_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N010-161128-4C-B-2X_S62_L001_RESULTS/N010-161128-4C-B-2X_S62_L001_Swarm_d1/N010-161128-4C-B-2X_S62_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N013-161205-5C-B_S52_L001_RESULTS/N013-161205-5C-B_S52_L001_Swarm_d1/N013-161205-5C-B_S52_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-1B-B_S56_L001_RESULTS/N09-161117-1B-B_S56_L001_Swarm_d1/N09-161117-1B-B_S56_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./Controlidentified.out/Controlidentified.out_Swarm_d1/Controlidentified.out_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N011-161130-2A-B_S74_L001_RESULTS/N011-161130-2A-B_S74_L001_Swarm_d1/N011-161130-2A-B_S74_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N013-161205-6A-B_S64_L001_RESULTS/N013-161205-6A-B_S64_L001_Swarm_d1/N013-161205-6A-B_S64_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-1B-F_S68_L001_RESULTS/N09-161117-1B-F_S68_L001_Swarm_d1/N09-161117-1B-F_S68_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./convert_barcodes.py/convert_barcodes.py_Swarm_d1/convert_barcodes.py_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N011-161130-2C-B_S86_L001_RESULTS/N011-161130-2C-B_S86_L001_Swarm_d1/N011-161130-2C-B_S86_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N02-160610-5D-F_S13_L001_RESULTS/N02-160610-5D-F_S13_L001_Swarm_d1/N02-160610-5D-F_S13_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-1C-B_S80_L001_RESULTS/N09-161117-1C-B_S80_L001_Swarm_d1/N09-161117-1C-B_S80_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./count_phytophthora_species.py/count_phytophthora_species.py_Swarm_d1/count_phytophthora_species.py_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N011-161130-3A-B_S3_L001_RESULTS/N011-161130-3A-B_S3_L001_Swarm_d1/N011-161130-3A-B_S3_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N02-160610-7C-F_S25_L001_RESULTS/N02-160610-7C-F_S25_L001_Swarm_d1/N02-160610-7C-F_S25_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-2A-B_S92_L001_RESULTS/N09-161117-2A-B_S92_L001_Swarm_d1/N09-161117-2A-B_S92_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./dada2_fastafiles/dada2_fastafiles_Swarm_d1/dada2_fastafiles_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N011-161130-3B-B_S15_L001_RESULTS/N011-161130-3B-B_S15_L001_Swarm_d1/N011-161130-3B-B_S15_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N03-160706-3C-F_S26_L001_RESULTS/N03-160706-3C-F_S26_L001_Swarm_d1/N03-160706-3C-F_S26_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-3A-F1-2x_S9_L001_RESULTS/N09-161117-3A-F1-2x_S9_L001_Swarm_d1/N09-161117-3A-F1-2x_S9_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GC1-0x_S82_L001_RESULTS/GC1-0x_S82_L001_Swarm_d1/GC1-0x_S82_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N011-161130-4AB_S88_L001_RESULTS/N011-161130-4AB_S88_L001_Swarm_d1/N011-161130-4AB_S88_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N03-160706-5C-F_S34_L001_RESULTS/N03-160706-5C-F_S34_L001_Swarm_d1/N03-160706-5C-F_S34_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-3B-B_S85_L001_RESULTS/N09-161117-3B-B_S85_L001_Swarm_d1/N09-161117-3B-B_S85_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GC1-100x_S11_L001_RESULTS/GC1-100x_S11_L001_Swarm_d1/GC1-100x_S11_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N01-160517-100-R_S17_L001_RESULTS/N01-160517-100-R_S17_L001_Swarm_d1/N01-160517-100-R_S17_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N03-160706-5D-F_S46_L001_RESULTS/N03-160706-5D-F_S46_L001_Swarm_d1/N03-160706-5D-F_S46_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-3B-F_S21_L001_RESULTS/N09-161117-3B-F_S21_L001_Swarm_d1/N09-161117-3B-F_S21_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GC1-10x_S94_L001_RESULTS/GC1-10x_S94_L001_Swarm_d1/GC1-10x_S94_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N01-160517-102-R_S65_L001_RESULTS/N01-160517-102-R_S65_L001_Swarm_d1/N01-160517-102-R_S65_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N03-160706-8A-F_S58_L001_RESULTS/N03-160706-8A-F_S58_L001_Swarm_d1/N03-160706-8A-F_S58_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-3C-B_S2_L001_RESULTS/N09-161117-3C-B_S2_L001_Swarm_d1/N09-161117-3C-B_S2_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GC2-0x_S23_L001_RESULTS/GC2-0x_S23_L001_Swarm_d1/GC2-0x_S23_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N01-160517-1PB-B_S5_L001_RESULTS/N01-160517-1PB-B_S5_L001_Swarm_d1/N01-160517-1PB-B_S5_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N04-160718-13B-F_S70_L001_RESULTS/N04-160718-13B-F_S70_L001_Swarm_d1/N04-160718-13B-F_S70_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-3C-F_S33_L001_RESULTS/N09-161117-3C-F_S33_L001_Swarm_d1/N09-161117-3C-F_S33_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GC2-100x_S47_L001_RESULTS/GC2-100x_S47_L001_Swarm_d1/GC2-100x_S47_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N01-160517-2WB-F_S1_L001_RESULTS/N01-160517-2WB-F_S1_L001_Swarm_d1/N01-160517-2WB-F_S1_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N05-161019-4B-B_S77_L001_RESULTS/N05-161019-4B-B_S77_L001_Swarm_d1/N05-161019-4B-B_S77_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-4A-F_S45_L001_RESULTS/N09-161117-4A-F_S45_L001_Swarm_d1/N09-161117-4A-F_S45_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GC2-10x_S35_L001_RESULTS/GC2-10x_S35_L001_Swarm_d1/GC2-10x_S35_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N01-160630-2B-F-2X_S29_L001_RESULTS/N01-160630-2B-F-2X_S29_L001_Swarm_d1/N01-160630-2B-F-2X_S29_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N05-161019-4C-F_S6_L001_RESULTS/N05-161019-4C-F_S6_L001_Swarm_d1/N05-161019-4C-F_S6_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-4B-B_S14_L001_RESULTS/N09-161117-4B-B_S14_L001_Swarm_d1/N09-161117-4B-B_S14_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GC3-0x_S59_L001_RESULTS/GC3-0x_S59_L001_Swarm_d1/GC3-0x_S59_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N01-160630-3A-F_S41_L001_RESULTS/N01-160630-3A-F_S41_L001_Swarm_d1/N01-160630-3A-F_S41_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N05-161019-5B-F_S37_L001_RESULTS/N05-161019-5B-F_S37_L001_Swarm_d1/N05-161019-5B-F_S37_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-4B-F_S57_L001_RESULTS/N09-161117-4B-F_S57_L001_Swarm_d1/N09-161117-4B-F_S57_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GC3-100x_S83_L001_RESULTS/GC3-100x_S83_L001_Swarm_d1/GC3-100x_S83_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N01-160630-41A-R_S89_L001_RESULTS/N01-160630-41A-R_S89_L001_Swarm_d1/N01-160630-41A-R_S89_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N06-161020-2C-F_S61_L001_RESULTS/N06-161020-2C-F_S61_L001_Swarm_d1/N06-161020-2C-F_S61_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-4C-B_S69_L001_RESULTS/N09-161117-4C-B_S69_L001_Swarm_d1/N09-161117-4C-B_S69_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GC3-10x_S71_L001_RESULTS/GC3-10x_S71_L001_Swarm_d1/GC3-10x_S71_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N01-160630-41D-R_S49_L001_RESULTS/N01-160630-41D-R_S49_L001_Swarm_d1/N01-160630-41D-R_S49_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N06-161020-5A-F_S30_L001_RESULTS/N06-161020-5A-F_S30_L001_Swarm_d1/N06-161020-5A-F_S30_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-4C-F_S81_L001_RESULTS/N09-161117-4C-F_S81_L001_Swarm_d1/N09-161117-4C-F_S81_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./generate_metapy.py/generate_metapy.py_Swarm_d1/generate_metapy.py_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N01-160630-42C-R_S20_L001_RESULTS/N01-160630-42C-R_S20_L001_Swarm_d1/N01-160630-42C-R_S20_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N06-161020-5B-B_S73_L001_RESULTS/N06-161020-5B-B_S73_L001_Swarm_d1/N06-161020-5B-B_S73_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-6A-F_S18_L001_RESULTS/N09-161117-6A-F_S18_L001_Swarm_d1/N09-161117-6A-F_S18_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GL4-0x_S95_L001_RESULTS/GL4-0x_S95_L001_Swarm_d1/GL4-0x_S95_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N01-160630-44A-R_S93_L001_RESULTS/N01-160630-44A-R_S93_L001_Swarm_d1/N01-160630-44A-R_S93_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N06-161020-5B-F_S42_L001_RESULTS/N06-161020-5B-F_S42_L001_Swarm_d1/N06-161020-5B-F_S42_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N09-161117-6B-F_S43_L001_RESULTS/N09-161117-6B-F_S43_L001_Swarm_d1/N09-161117-6B-F_S43_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GL4-100x_S24_L001_RESULTS/GL4-100x_S24_L001_Swarm_d1/GL4-100x_S24_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N01-160630-45A-R_S10_L001_RESULTS/N01-160630-45A-R_S10_L001_Swarm_d1/N01-160630-45A-R_S10_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N06-161020-5C-F_S54_L001_RESULTS/N06-161020-5C-F_S54_L001_Swarm_d1/N06-161020-5C-F_S54_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./Rplots.pdf/Rplots.pdf_Swarm_d1/Rplots.pdf_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GL4-10x_S12_L001_RESULTS/GL4-10x_S12_L001_Swarm_d1/GL4-10x_S12_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N01-160630-45C-R_S22_L001_RESULTS/N01-160630-45C-R_S22_L001_Swarm_d1/N01-160630-45C-R_S22_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N06-161020-6A-B_S66_L001_RESULTS/N06-161020-6A-B_S66_L001_Swarm_d1/N06-161020-6A-B_S66_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./swarm_results/swarm_results_Swarm_d1/swarm_results_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GL5-0x_S36_L001_RESULTS/GL5-0x_S36_L001_Swarm_d1/GL5-0x_S36_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N01-160630-6A-F-2X_S53_L001_RESULTS/N01-160630-6A-F-2X_S53_L001_Swarm_d1/N01-160630-6A-F-2X_S53_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N06-161020-6A-F_S78_L001_RESULTS/N06-161020-6A-F_S78_L001_Swarm_d1/N06-161020-6A-F_S78_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./Undetermined_S0_L001_RESULTS/Undetermined_S0_L001_Swarm_d1/Undetermined_S0_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GL5-100x_S60_L001_RESULTS/GL5-100x_S60_L001_Swarm_d1/GL5-100x_S60_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N012-161201-2A-B_S27_L001_RESULTS/N012-161201-2A-B_S27_L001_Swarm_d1/N012-161201-2A-B_S27_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N06-161020-6B-F_S90_L001_RESULTS/N06-161020-6B-F_S90_L001_Swarm_d1/N06-161020-6B-F_S90_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./GL5-10x_S48_L001_RESULTS/GL5-10x_S48_L001_Swarm_d1/GL5-10x_S48_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N012-161201-2B-B_S39_L001_RESULTS/N012-161201-2B-B_S39_L001_Swarm_d1/N012-161201-2B-B_S39_L001_swarm_results_1.RESULTS ./swarm_results/
 
cp ./N06-161020-7A-F_S7_L001_RESULTS/N06-161020-7A-F_S7_L001_Swarm_d1/N06-161020-7A-F_S7_L001_swarm_results_1.RESULTS ./swarm_results/




################################################################################
# bowtie:


cp ./1_Phytophthora_cactorum_AF266772identified.out/1_Phytophthora_cactorum_AF266772identified.out_bowtie/1_Phytophthora_cactorum_AF266772identified.out_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GL6-0x_S72_L001_RESULTS/GL6-0x_S72_L001_bowtie/GL6-0x_S72_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N012-161201-5A-B_S51_L001_RESULTS/N012-161201-5A-B_S51_L001_bowtie/N012-161201-5A-B_S51_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N06-161020-7B-F_S19_L001_RESULTS/N06-161020-7B-F_S19_L001_bowtie/N06-161020-7B-F_S19_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./8c_Phytophthora_lateralis_AF266804identified.out/8c_Phytophthora_lateralis_AF266804identified.out_bowtie/8c_Phytophthora_lateralis_AF266804identified.out_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GL6-100x_S96_L001_RESULTS/GL6-100x_S96_L001_bowtie/GL6-100x_S96_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N012-161201-5B-B_S63_L001_RESULTS/N012-161201-5B-B_S63_L001_bowtie/N012-161201-5B-B_S63_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N07-161028-1B-F_S31_L001_RESULTS/N07-161028-1B-F_S31_L001_bowtie/N07-161028-1B-F_S31_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./8c_Phytophthora_ramorum_JQ653034identified.out/8c_Phytophthora_ramorum_JQ653034identified.out_bowtie/8c_Phytophthora_ramorum_JQ653034identified.out_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GL6-10x_S84_L001_RESULTS/GL6-10x_S84_L001_bowtie/GL6-10x_S84_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N012-161201-6B-B_S75_L001_RESULTS/N012-161201-6B-B_S75_L001_bowtie/N012-161201-6B-B_S75_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N07-161028-3B-F_S55_L001_RESULTS/N07-161028-3B-F_S55_L001_bowtie/N07-161028-3B-F_S55_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./8d_Phytophthora_austrocedri_DQ995184identified.out/8d_Phytophthora_austrocedri_DQ995184identified.out_bowtie/8d_Phytophthora_austrocedri_DQ995184identified.out_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./interogate_results.py/interogate_results.py_bowtie/interogate_results.py_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N012-161201-6C-B_S87_L001_RESULTS/N012-161201-6C-B_S87_L001_bowtie/N012-161201-6C-B_S87_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N07-161028-3C-B_S67_L001_RESULTS/N07-161028-3C-B_S67_L001_bowtie/N07-161028-3C-B_S67_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./all_results_files.txt/all_results_files.txt_bowtie/all_results_files.txt_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./log_files/log_files_bowtie/log_files_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N013-161205-10B-B_S76_L001_RESULTS/N013-161205-10B-B_S76_L001_bowtie/N013-161205-10B-B_S76_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N07-161028-3C-F_S79_L001_RESULTS/N07-161028-3C-F_S79_L001_bowtie/N07-161028-3C-F_S79_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./bowtie_results/bowtie_results_bowtie/bowtie_results_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./metapy.sh/metapy.sh_bowtie/metapy.sh_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N013-161205-1C-B_S4_L001_RESULTS/N013-161205-1C-B_S4_L001_bowtie/N013-161205-1C-B_S4_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N07-161028-5A-F_S91_L001_RESULTS/N07-161028-5A-F_S91_L001_bowtie/N07-161028-5A-F_S91_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./Control_C1identified.out/Control_C1identified.out_bowtie/Control_C1identified.out_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./mv.sh/mv.sh_bowtie/mv.sh_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N013-161205-3A-B_S16_L001_RESULTS/N013-161205-3A-B_S16_L001_bowtie/N013-161205-3A-B_S16_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N07-161028-5C-F_S8_L001_RESULTS/N07-161028-5C-F_S8_L001_bowtie/N07-161028-5C-F_S8_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./Control_C2identified.out/Control_C2identified.out_bowtie/Control_C2identified.out_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N010-161128-2A-F_S38_L001_RESULTS/N010-161128-2A-F_S38_L001_bowtie/N010-161128-2A-F_S38_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N013-161205-5A-B_S28_L001_RESULTS/N013-161205-5A-B_S28_L001_bowtie/N013-161205-5A-B_S28_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-1A-B_S32_L001_RESULTS/N09-161117-1A-B_S32_L001_bowtie/N09-161117-1A-B_S32_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./Control_C3identified.out/Control_C3identified.out_bowtie/Control_C3identified.out_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N010-161128-4A-F_S50_L001_RESULTS/N010-161128-4A-F_S50_L001_bowtie/N010-161128-4A-F_S50_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N013-161205-5B-B_S40_L001_RESULTS/N013-161205-5B-B_S40_L001_bowtie/N013-161205-5B-B_S40_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-1A-F_S44_L001_RESULTS/N09-161117-1A-F_S44_L001_bowtie/N09-161117-1A-F_S44_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./Control_C4identified.out/Control_C4identified.out_bowtie/Control_C4identified.out_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N010-161128-4C-B-2X_S62_L001_RESULTS/N010-161128-4C-B-2X_S62_L001_bowtie/N010-161128-4C-B-2X_S62_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N013-161205-5C-B_S52_L001_RESULTS/N013-161205-5C-B_S52_L001_bowtie/N013-161205-5C-B_S52_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-1B-B_S56_L001_RESULTS/N09-161117-1B-B_S56_L001_bowtie/N09-161117-1B-B_S56_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./Controlidentified.out/Controlidentified.out_bowtie/Controlidentified.out_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N011-161130-2A-B_S74_L001_RESULTS/N011-161130-2A-B_S74_L001_bowtie/N011-161130-2A-B_S74_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N013-161205-6A-B_S64_L001_RESULTS/N013-161205-6A-B_S64_L001_bowtie/N013-161205-6A-B_S64_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-1B-F_S68_L001_RESULTS/N09-161117-1B-F_S68_L001_bowtie/N09-161117-1B-F_S68_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./convert_barcodes.py/convert_barcodes.py_bowtie/convert_barcodes.py_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N011-161130-2C-B_S86_L001_RESULTS/N011-161130-2C-B_S86_L001_bowtie/N011-161130-2C-B_S86_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N02-160610-5D-F_S13_L001_RESULTS/N02-160610-5D-F_S13_L001_bowtie/N02-160610-5D-F_S13_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-1C-B_S80_L001_RESULTS/N09-161117-1C-B_S80_L001_bowtie/N09-161117-1C-B_S80_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./count_phytophthora_species.py/count_phytophthora_species.py_bowtie/count_phytophthora_species.py_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N011-161130-3A-B_S3_L001_RESULTS/N011-161130-3A-B_S3_L001_bowtie/N011-161130-3A-B_S3_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N02-160610-7C-F_S25_L001_RESULTS/N02-160610-7C-F_S25_L001_bowtie/N02-160610-7C-F_S25_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-2A-B_S92_L001_RESULTS/N09-161117-2A-B_S92_L001_bowtie/N09-161117-2A-B_S92_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./dada2_fastafiles/dada2_fastafiles_bowtie/dada2_fastafiles_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N011-161130-3B-B_S15_L001_RESULTS/N011-161130-3B-B_S15_L001_bowtie/N011-161130-3B-B_S15_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N03-160706-3C-F_S26_L001_RESULTS/N03-160706-3C-F_S26_L001_bowtie/N03-160706-3C-F_S26_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-3A-F1-2x_S9_L001_RESULTS/N09-161117-3A-F1-2x_S9_L001_bowtie/N09-161117-3A-F1-2x_S9_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GC1-0x_S82_L001_RESULTS/GC1-0x_S82_L001_bowtie/GC1-0x_S82_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N011-161130-4AB_S88_L001_RESULTS/N011-161130-4AB_S88_L001_bowtie/N011-161130-4AB_S88_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N03-160706-5C-F_S34_L001_RESULTS/N03-160706-5C-F_S34_L001_bowtie/N03-160706-5C-F_S34_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-3B-B_S85_L001_RESULTS/N09-161117-3B-B_S85_L001_bowtie/N09-161117-3B-B_S85_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GC1-100x_S11_L001_RESULTS/GC1-100x_S11_L001_bowtie/GC1-100x_S11_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N01-160517-100-R_S17_L001_RESULTS/N01-160517-100-R_S17_L001_bowtie/N01-160517-100-R_S17_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N03-160706-5D-F_S46_L001_RESULTS/N03-160706-5D-F_S46_L001_bowtie/N03-160706-5D-F_S46_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-3B-F_S21_L001_RESULTS/N09-161117-3B-F_S21_L001_bowtie/N09-161117-3B-F_S21_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GC1-10x_S94_L001_RESULTS/GC1-10x_S94_L001_bowtie/GC1-10x_S94_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N01-160517-102-R_S65_L001_RESULTS/N01-160517-102-R_S65_L001_bowtie/N01-160517-102-R_S65_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N03-160706-8A-F_S58_L001_RESULTS/N03-160706-8A-F_S58_L001_bowtie/N03-160706-8A-F_S58_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-3C-B_S2_L001_RESULTS/N09-161117-3C-B_S2_L001_bowtie/N09-161117-3C-B_S2_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GC2-0x_S23_L001_RESULTS/GC2-0x_S23_L001_bowtie/GC2-0x_S23_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N01-160517-1PB-B_S5_L001_RESULTS/N01-160517-1PB-B_S5_L001_bowtie/N01-160517-1PB-B_S5_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N04-160718-13B-F_S70_L001_RESULTS/N04-160718-13B-F_S70_L001_bowtie/N04-160718-13B-F_S70_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-3C-F_S33_L001_RESULTS/N09-161117-3C-F_S33_L001_bowtie/N09-161117-3C-F_S33_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GC2-100x_S47_L001_RESULTS/GC2-100x_S47_L001_bowtie/GC2-100x_S47_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N01-160517-2WB-F_S1_L001_RESULTS/N01-160517-2WB-F_S1_L001_bowtie/N01-160517-2WB-F_S1_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N05-161019-4B-B_S77_L001_RESULTS/N05-161019-4B-B_S77_L001_bowtie/N05-161019-4B-B_S77_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-4A-F_S45_L001_RESULTS/N09-161117-4A-F_S45_L001_bowtie/N09-161117-4A-F_S45_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GC2-10x_S35_L001_RESULTS/GC2-10x_S35_L001_bowtie/GC2-10x_S35_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N01-160630-2B-F-2X_S29_L001_RESULTS/N01-160630-2B-F-2X_S29_L001_bowtie/N01-160630-2B-F-2X_S29_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N05-161019-4C-F_S6_L001_RESULTS/N05-161019-4C-F_S6_L001_bowtie/N05-161019-4C-F_S6_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-4B-B_S14_L001_RESULTS/N09-161117-4B-B_S14_L001_bowtie/N09-161117-4B-B_S14_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GC3-0x_S59_L001_RESULTS/GC3-0x_S59_L001_bowtie/GC3-0x_S59_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N01-160630-3A-F_S41_L001_RESULTS/N01-160630-3A-F_S41_L001_bowtie/N01-160630-3A-F_S41_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N05-161019-5B-F_S37_L001_RESULTS/N05-161019-5B-F_S37_L001_bowtie/N05-161019-5B-F_S37_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-4B-F_S57_L001_RESULTS/N09-161117-4B-F_S57_L001_bowtie/N09-161117-4B-F_S57_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GC3-100x_S83_L001_RESULTS/GC3-100x_S83_L001_bowtie/GC3-100x_S83_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N01-160630-41A-R_S89_L001_RESULTS/N01-160630-41A-R_S89_L001_bowtie/N01-160630-41A-R_S89_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N06-161020-2C-F_S61_L001_RESULTS/N06-161020-2C-F_S61_L001_bowtie/N06-161020-2C-F_S61_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-4C-B_S69_L001_RESULTS/N09-161117-4C-B_S69_L001_bowtie/N09-161117-4C-B_S69_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GC3-10x_S71_L001_RESULTS/GC3-10x_S71_L001_bowtie/GC3-10x_S71_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N01-160630-41D-R_S49_L001_RESULTS/N01-160630-41D-R_S49_L001_bowtie/N01-160630-41D-R_S49_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N06-161020-5A-F_S30_L001_RESULTS/N06-161020-5A-F_S30_L001_bowtie/N06-161020-5A-F_S30_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-4C-F_S81_L001_RESULTS/N09-161117-4C-F_S81_L001_bowtie/N09-161117-4C-F_S81_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./generate_metapy.py/generate_metapy.py_bowtie/generate_metapy.py_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N01-160630-42C-R_S20_L001_RESULTS/N01-160630-42C-R_S20_L001_bowtie/N01-160630-42C-R_S20_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N06-161020-5B-B_S73_L001_RESULTS/N06-161020-5B-B_S73_L001_bowtie/N06-161020-5B-B_S73_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-6A-F_S18_L001_RESULTS/N09-161117-6A-F_S18_L001_bowtie/N09-161117-6A-F_S18_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GL4-0x_S95_L001_RESULTS/GL4-0x_S95_L001_bowtie/GL4-0x_S95_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N01-160630-44A-R_S93_L001_RESULTS/N01-160630-44A-R_S93_L001_bowtie/N01-160630-44A-R_S93_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N06-161020-5B-F_S42_L001_RESULTS/N06-161020-5B-F_S42_L001_bowtie/N06-161020-5B-F_S42_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N09-161117-6B-F_S43_L001_RESULTS/N09-161117-6B-F_S43_L001_bowtie/N09-161117-6B-F_S43_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GL4-100x_S24_L001_RESULTS/GL4-100x_S24_L001_bowtie/GL4-100x_S24_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N01-160630-45A-R_S10_L001_RESULTS/N01-160630-45A-R_S10_L001_bowtie/N01-160630-45A-R_S10_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N06-161020-5C-F_S54_L001_RESULTS/N06-161020-5C-F_S54_L001_bowtie/N06-161020-5C-F_S54_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./Rplots.pdf/Rplots.pdf_bowtie/Rplots.pdf_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GL4-10x_S12_L001_RESULTS/GL4-10x_S12_L001_bowtie/GL4-10x_S12_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N01-160630-45C-R_S22_L001_RESULTS/N01-160630-45C-R_S22_L001_bowtie/N01-160630-45C-R_S22_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N06-161020-6A-B_S66_L001_RESULTS/N06-161020-6A-B_S66_L001_bowtie/N06-161020-6A-B_S66_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./bowtie_results/swarm_results_bowtie/swarm_results_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GL5-0x_S36_L001_RESULTS/GL5-0x_S36_L001_bowtie/GL5-0x_S36_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N01-160630-6A-F-2X_S53_L001_RESULTS/N01-160630-6A-F-2X_S53_L001_bowtie/N01-160630-6A-F-2X_S53_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N06-161020-6A-F_S78_L001_RESULTS/N06-161020-6A-F_S78_L001_bowtie/N06-161020-6A-F_S78_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./Undetermined_S0_L001_RESULTS/Undetermined_S0_L001_bowtie/Undetermined_S0_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GL5-100x_S60_L001_RESULTS/GL5-100x_S60_L001_bowtie/GL5-100x_S60_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N012-161201-2A-B_S27_L001_RESULTS/N012-161201-2A-B_S27_L001_bowtie/N012-161201-2A-B_S27_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N06-161020-6B-F_S90_L001_RESULTS/N06-161020-6B-F_S90_L001_bowtie/N06-161020-6B-F_S90_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./GL5-10x_S48_L001_RESULTS/GL5-10x_S48_L001_bowtie/GL5-10x_S48_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N012-161201-2B-B_S39_L001_RESULTS/N012-161201-2B-B_S39_L001_bowtie/N012-161201-2B-B_S39_L001_bowtie_perfect.RESULTS ./bowtie_results/
 
cp ./N06-161020-7A-F_S7_L001_RESULTS/N06-161020-7A-F_S7_L001_bowtie/N06-161020-7A-F_S7_L001_bowtie_perfect.RESULTS ./bowtie_results/

########################################
# summerize:
cd bowtie_results
grep -H "percent of assembled-reads clustering with database" * > bowtie_percentage_cluster_with_db.txt
grep -H "Total number of assembled sequences" * > total_number_of_assembled_sequences.txt

cd ../
cd swarm_results
grep -H "percent of assembled-reads clustering with database" * > swarm_percentage_cluster_with_db.txt
grep -H "Total number of assembled sequences" * > total_number_of_assembled_sequences.txt


