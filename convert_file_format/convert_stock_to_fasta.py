from Bio import SeqIO
records = SeqIO.parse("G_B-M_backontospry.sto", "stockholm")
count = SeqIO.write(records, "Gpal_m_B_back_onto_spry_definition_hmmre_reformatted.fasta", "fasta")
print "Converted %i records" % count
