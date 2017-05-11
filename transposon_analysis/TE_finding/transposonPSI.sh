cd '/home/peter/.gvfs/SFTP for pt40963 on 143.234.96.240/mnt/shared/users/pt40963/nematode/Mitochondria_20160311/psi_transposon_finding' 

fasta_filenames=../*.fa

for f in ${fasta_filenames}
do
	echo "Running  ${f}"
	#step1
	cmd="perl /home/peter/Downloads/TransposonPSI_08222010/transposonPSI.pl ${f} nuc" 
	echo ${cmd}
	eval ${cmd}
	wait
done




