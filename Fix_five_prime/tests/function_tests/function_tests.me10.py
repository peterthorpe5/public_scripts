def complement(DNA):
    DNA_1 = DNA.replace("A", "t")
    DNA_2 = DNA_1.replace("T", "a")
    DNA_3 = DNA_2.replace("C", "g")
    complement = DNA_3.replace("G", "c")
    #print complement
    return complement.upper()

def reverse_complement(DNA):
    compl = complement(DNA)
    return compl[::-1]


def find_positive_next_ATG(transcriptome_record, position, strand):
    "function to find the next ATG"
    #print ("start position = %d" % position)
    transcriptome_record= transcriptome_record[position:]
    if strand =="+":
        next_codon = ""
        start = 0
        end = 3
        for i in range(len(transcriptome_record)-60):
            start += 3
            end += 3
            next_codon = transcriptome_record[start:end]
            print("start %i, end %i, codon: %s" % (start, end, next_codon))
            #print next_codon
            #next_codon = transcript.seq[position+3:position+6]
            if next_codon == "ATG":
                #print "fucking yes"
                next_methionine = start
                return next_methionine + position  # Note restoring offet!
        #return None
        raise ValueError("No start codon found")
    else:
        raise ValueError("Called on strand %r" % strand)


def find_downstream_start(transcript, current_start, strand):
    if strand == "+":
        print("Looking for ATG after %d in sequence %s" % (current_start, transcript))
        if transcript[current_start:current_start+3] != "ATG":
            print("WARNING - existing annotation does not start ATG")
        return find_positive_next_ATG(transcript, current_start, strand)
    elif strand == "-":
        print("Looking for CAT (i.e. ATG rev-comp) before %d in sequence %s" % (current_start, transcript))
        new = find_positive_next_ATG(reverse_complement(transcript),
                                     len(transcript) - current_start, "+")
        assert new is not None
        return len(transcript) - new
    else:
        raise ValueError("Bad strand value %r" % strand)
 
#print reverse_complement("ATG")

me10_tran = """GGTTACATTTGATTTATAGAAAAAAATAACATACATCATAGTCAAAAACAAAACATAGAA
TTACAATAAGTTTTATGCTACAACAGGTTGGGCTACGGTAGGGGACATTGATTTAAGGAT
ATTTTCATGCATATCATCTAGAGCATCTTTATAATTGTCCAAGTCATCCATTTGATGGTA
AACTTTTTGTGGATCTGAGGTTTGATAGTTGATTGTGGGCTCAAATTCATTTGTTGCCTT
GAGATATTGTATAACTTTTTGGAACCACTTGTTTTTGGACTTATTAAGTTCATTCCTTGC
TCTTTCCACAGCAGAATAGTCAACTTTTCCAGCAAAATATTCTCTTTGCATGTCGTTTAA
TGTCTGCGAATCGTCCATCAAGTTCATCTCAGCCTCTCCCAATTCGTACAAACTAGCTTT
TACTGCACAGAAATCTTGGTCTATTAATGGTTGCATCATTAATGCATTGCTTTCTAGTAC
ATTGTAGCTTATGATGACCACCACTAAACTGGAAAACATTGCGATTTCCTTGAAAGCCAT
ATTTGGTCGGAGTTTCCGCAAGTGTATTAGATATTTA""".replace("\n", "")


me10 = """AAATATCTAATACACTTGCGGAAACTCCGACCAAATATGGCTTTCAAGGAAATCGCAATG
TTTTCCAGTTTAGTGGTGGTCATCATAAGCTACAATGTACTAGAAAGCAATGCATTAATG
ATGCAACCATTAATAGACCAAGATTTCTGTGCAGTAAAAGCTAGTTTGTACGAATTGGGA
GAGGCTGAGATGAACTTGATGGACGATTCGCAGACATTAAACGACATGCAAAGAGAATAT
TTTGCTGGAAAAGTTGACTATTCTGCTGTGGAAAGAGCAAGGAATGAACTTAATAAGTCC
AAAAACAAGTGGTTCCAAAAAGTTATACAATATCTCAAGGCAACAAATGAATTTGAGCCC
ACAATCAACTATCAAACCTCAGATCCACAAAAAGTTTACCATCAAATGGATGACTTGGAC
AATTATAAAGATGCTCTAGATGATATGCATGAAAATATCCTTAAATCAATGTCCCCTACC
GTAGCCCAACCTGTTGTAGCATAA""".replace("\n", "")

#print transcript[140:]

#transcript = "ATGzzzGGGTTTbbbATGeeeCC"

#broken_up_for_reading = "ATG NNN GGG TTT NNN ATG GGCCC"


me10_rc = reverse_complement(me10) # CDS


print "Doing it on the reversed strand so its all positive (easy case)"
start = reverse_complement(me10_tran).index(me10)
print reverse_complement(me10_tran)
print ("-" * start) + me10
print find_downstream_start(reverse_complement(me10_tran), start, "+")

print "Doing it on the reverse strand as is (hard case)"
print find_downstream_start(me10_tran, len(me10_tran) - reverse_complement(me10_tran).index(me10), "-")
print "..."

end_pos_negative = me10_tran.find(me10_rc)
end_pos_negative = me10_tran.index(me10_rc)


start_pos_negative = (end_pos_negative+len(me10))

print "start_pos_negative = ", start_pos_negative, "end_pos_negative = ", end_pos_negative

#start = find_negative_next_ATG(me10_tran,end_pos_negative, start_pos_negative,"-")
start = find_downstream_start(me10_tran,end_pos_negative, start_pos_negative,"-")


print "start = ", start

print "NEW gene"
print reverse_complement(me10_tran[end_pos_negative:start])

print "next start = ", start

#print me10_tran[start_pos:end]
