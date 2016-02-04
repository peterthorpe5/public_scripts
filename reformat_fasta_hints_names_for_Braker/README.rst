READme info for script reformat fasta and hints names
======================================================

basic usage:

``./re_format_fasta_hints.py`` -h 


what?
script to reformat fasta and hints file names.

why?
Braker (GeneMark) fails with either "|" or long names. e.g. scaffold13578|length67900
even though Braker will rename these scaffolds for GeneMark - it will fail with 
GeneMark version 4.29.
(Augustus3.1, Bamtools 2.4, GeneMark-ES / ET v.4.32  ( the Readme says to install 4.32, but I cannot find this version online - GeneMark-ET 4.29)

how?
scaffold13578|length67900 will be replace with --prefix  e.g Rp

scaffold13578|length67900  >   Rp67900    = this will now work. 

author: Peter Thorpe September 2015. The James Hutton Insitute, Dundee, UK.

example Braker commands (generate hints):
 ``bam2hints --intronsonly --in=rnaseq.fs.bam --out=hints.gff``

reformat names:

 `` python ~/misc_python/re_format_fasta_hints.py --hints ./all_RNAseq_mapped/hints.gff -f Myzus_cerasi_genome_assembly.v1.fasta --prefix Mc -o Mc_v1 ``
 
 run Braker
 
 `` braker.pl --hints=Mc_v1_alt.hints --AUGUSTUS_CONFIG_PATH=/home/peter/Downloads/augustus-3.1/bin --species=pea_aphid --BAMTOOLS_PATH=/home/peter/Downloads/bamtools-master/bin --cores 1 --optCfgFile=extrinsic.bug.cfg --UTR on --gff3 --GENEMARK_PATH=/home/peter/Downloads/gm_et_linux_64/gmes_petap --workingdir=/home/peter/Desktop/Mc_braker --genome=Mc_v1_alt.fasta ``



 
 
 
 
OTHER ERRORS ENCOUNTED AND FIXED BY:
====================================

Generate  hints file by using:
 ``bam2hints --intronsonly --in=rnaseq.fs.bam --out=hints.gff``
 
This input failed as did the bam as input to Braker, which all failed….. So I did some investigating. I think I have found a solution which may help your future users. Braker has been running error free for hours now (it is still not finished – but has not crashed):
 
As everything failed with GeneMark, I tried running this program alone:
 
gmes_petap.pl --ET hints.gff --sequence R.padi_final_genome.v1.fasta --et_score 10
 
The error came back as:
 
"inappropriate ioctl for device" …
Google serach:
(http://stackoverflow.com/questions/1605195/inappropriate-ioctl-for-device)

which identified this being this problem
GLIBCXX_3.4.20

(http://www.unix.com/shell-programming-and-scripting/20812-inappropriate-ioctl-device.html)

How to fix:

http://askubuntu.com/questions/575505/glibcxx-3-4-20-not-found-how-to-fix-this-error

sudo apt-get install libstdc++6

This didnt update anything as the latest was already intalled. 

sudo add-apt-repository ppa:ubuntu-toolchain-r/test 
sudo apt-get update
sudo apt-get upgrade
sudo apt-get dist-upgrade

these fixed the problem.
 
I hope this is helpful to other people
 
Cheers,
 
Peter Thorpe
