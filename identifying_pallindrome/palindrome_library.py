"""Library functions for palidrome (simple inverted repeats, SIR) searching.

version - changes
=================
001 - 10th April - this code is optimised for speed and not clarity!!!

Pete the python tamer optimised this!

Function to find all SIR of a given length in a larger sequence. with
maximium probability of being random
"""
from sys import stdin,argv
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
import time
from string import maketrans
import random

#Ignores ambiguous letters like D
complement_trans_table = maketrans("ACGTacgt", "TGCAtgca") 

def load_genome_sequence(filename):
    """this is a function to open a desired file and return the contents as a
    string - for speed purposes"""
    file_handle = open(filename, "r")
    record = SeqIO.read(file_handle, "fasta")
    file_handle.close()
    return str(record.seq)

def write_shuffled_genome(original_filename, shuffled_filename) :
    """Write a shuffled version of the sequence to a FASTA file."""

    record = SeqIO.read(open(original_filename), "fasta")

    data = list(record.seq)
    random.shuffle(data)
    shuffled_string = "".join(data)

    f= open(shuffled_filename, 'w')
    print >> f, ">random version of %s" % record.id
    #Elegant approach to line wrapping, writing out line by line@
    for index in range(0, len(shuffled_string), 70) :
        print >> f, shuffled_string[index:index+70]
    f.close()

def count_mismatches(seq, half_len):
    """ When passed a str(Bio.Seq seq), representing a nucleotide sequence,
        treats it as a short inverted repeat, and returns the number of 
        mismatched compared to its reverse complement for half the length
        of the sequence
    """    
    global complement_trans_table
    #assert half_len == len(seq)/2 #TODO - for testing, this is SLOW!
    #half_len = len(seq)/2
    #Want to use the python string's translate method to map A->T etc.
    second_half_seq_comp = seq[-half_len:].translate(complement_trans_table)
    mismatches = 0
    for i in range(half_len) :
        if seq[i] != second_half_seq_comp[-i-1] :
            mismatches += 1
    return mismatches

assert 0 == count_mismatches("AAAAATTTTT", 5)
assert 5 == count_mismatches("AAAAACCCCC", 5)


def base_composition(sequence):
    """function to work out the 'C,G,A,T' content of the genome/sequecne
    in question. The function takes a string or Bio.Seq Seq as the input"""
    s1 = str(sequence).upper()
    length = float(len(s1))
    bc= {'a': float(s1.count('A')/ length),
         't': float(s1.count('T')/ length),
         'g': float(s1.count('G')/ length),
         'c': float(s1.count('C')/ length)}
    return bc

def comb(N,k):
    """Combinations of N things taken k at a time.
    Exact long integer is computed.  Based on function from scipy.

    Notes:
      - If k > N, N < 0, or k < 0, then a 0 is returned.
    """
    if (k > N) or (N < 0) or (k < 0):
        return 0L
    val = 1L
    for j in xrange(min(k, N-k)):
        val = (val*(N-j))//(j+1)
    return val

def prob(full_length, mismatches, p_match) :
    """ this function is used to work out the probability of a palindrome of
    length - 'length', with x amount of mismatches, based on the probability of
    getting a match (p_match)

    the combination funtctio (comb) is used instead of scipy to work out the number
    of possible ways of getting that mismatch and length combination to work out
    the probability of that occuring at random

    probability is defined by:
    p = combinatons of mismatche (N choose K) * (probability of a match based on the
    codon bias ** number of matches) * (probability of a mismatch ** number of
    mismatches)

    (** to the power of)
    """
    halflength= full_length //2
    matches= halflength-mismatches
    combinations_of_mismatch = comb(halflength,mismatches)
    p_mismatch= 1-p_match
    return combinations_of_mismatch*\
           ((p_match**matches)*(p_mismatch**mismatches))

def base_composition_sequence_probability(dna, bc):
    """this is a function to return the probability of a given DNA sequence
    occuring at random. The base composition of the given string is taken
    into consideration(bc). dna should be a Bio.Seq seq. The given DNA sequence is
    treated as a pallindrome
    for a given sequence we need to find if its opposite palindromic pair
    is a match of miss match"""
    #we only need to use half of the length of the sequence so...
    halflength= (len(dna)/2)
    #call in function to return mismatches in the sequence
    mismatches=count_mismatches(dna, halflength)
    #work out the number of matches in a sequence
    matches= halflength-mismatches
    #we need to make a list of all the matching and mismatching pairs
    #the number of ways of getting that num of mismatches in a given sequence
    combinations_of_mismatch = comb(halflength,mismatches)
    #update values of bc using function
    #bc=base_composition(dna)
    #print bc
    matching_bases= ['at', 'ta', 'gc', 'cg']
    p_match = sum([bc[m[0]] * bc[m[1]] for m in matching_bases])
    p_mismatch= 1-p_match
    #the formula to work out the probability of a given sequence occuring
    #randomly in a genome.
    probability= combinations_of_mismatch*\
                 ((p_match**matches)*(p_mismatch**mismatches))
    #print 'pmatch', p_match
    return probability, mismatches
                                   
def iteration_length(data, lower, higher, probability_threshold=None):
    """this is a funtion to define the length of SIR we are looking for,

    lower to higher is the range of lengths we are interested in and should be an
    int.

    data is a string or Bio.Seq sequence.
    
    function should return sequences that only fullfill the limits defined by
    is_palindrome function, the lenth of this seq and the start position.
    """
    global complement_trans_table
    data_comp = data.translate(complement_trans_table)

    #call function to work out the codon bias to work out the probability
    #based on the codon bias of the whole genome...
    bc = base_composition(data)

    #define maches bases
    matching_bases= ['at', 'ta', 'gc', 'cg']
    #probability of a match... work out based on codon bias (bc)
    p_match = 0
    for i in matching_bases:
        p_match += bc[i[0]] * bc[i[1]]

    assert p_match == sum([bc[m[0]] * bc[m[1]] for m in matching_bases])

    #Big speed boost, cache the probabilities - these don't need to look at
    #actual potential palindrome sequence.
    #create dictionary of possible probability using length and mismatches
    #as the indexes - look up mismatches and length, returns the probability
    prob_dict = {}
    for length in range(lower, higher):
        for mismatches in range(0, (length//2)+1) :
            #call prob function to work out the probability... put in dict
            prob_dict[(length, mismatches)] = prob(length, mismatches, p_match)

    #create an empty list to put the result in
    palindrome_list = []
    for length in range(lower, higher):
        #print "Searching for length %i" % length
        #defined mismatch threshold - - - - IMPORTANT!!!!
        half_length = length // 2
        mismatches_threshold = 0.5*half_length
        half_length_list = range(half_length)
        for start in range(0, (len(data)-length+1)):
            seq = data[start:start+length]
            seq_comp = data_comp[start:start+length]
            mismatches = 0
            for i in half_length_list :
                if seq[i] != seq_comp[-i-1] :
                    mismatches += 1
            if mismatches <= mismatches_threshold :
                #Only now look up the probability,
                probability = prob_dict[(length, mismatches)]

                #assert probability, missmatches == base_composition_sequence_probability(seq, bc)

                if probability <= probability_threshold :
                    #results into list, as a list of tuples
                    palindrome_list+=[(length,start,probability,mismatches,seq)]
    return palindrome_list

def write_palindromes_to_file(outfilename,output):
    """ this is a function to write the results of a program to a defined file"""
    #this opens the destined file to write 'what we want' to it ('w')
    fh= open(outfilename, "w")
    #the output needs to be formatted in a way that we want it
    #this print statement write to out outfile

    #for loop to convert the individual tuples(i) into formated string.
    #The formated strings need to be added to the final list...
    #then printed to the desired file.
    print >> fh,"#length\tposition\tprobability\tmismatches\tpalindrome_seq"
    for i in output:
        output_formated= '%d\t%d\t%.3e\t%d\t%s' %(i[0], i[1],i[2],i[3],i[4])
        print >> fh, output_formated
        #this print statement write to out outfile
    
    #closes the outfile
    fh.close()
    return None #replace with None when happy with output!
    
########################################################################
### Testing!
########################################################################

#Example 1
#=========
seq = "CCCCCCCAGGGGGGG"
probability, mismatches = base_composition_sequence_probability(seq, base_composition(seq))
assert mismatches == 0
assert probability == 0.0029737530677523246
#print bc
palindrome_list = iteration_length(seq, 15,150, 1.0/len(seq)) #,fh
assert palindrome_list == [(len(seq), 0, probability, mismatches, seq)]
del seq, probability, mismatches, palindrome_list


#Example 1a
#=========
seq1 = "TAAAACCCCCTGGGGGTTTT"
bc = base_composition(seq1)
probability, mismatches = base_composition_sequence_probability(seq1, base_composition(seq1))
assert mismatches == 3
#assert probability == 0.0029737530677523246
palindrome_list = iteration_length(seq1, 19,20, 1.0/len(seq1)) #,fh
#assert palindrome_list == [(len(seq), 1, probability, mismatches, seq)]
del seq1, bc, probability, mismatches, palindrome_list

###Example 2
#=========
#going to search using lengths 15 upwards...
seq = "CCCCCCCAGGGGGGGG"
palindrome_list = []

#Expect this palindrome at position 0,
pal = "CCCCCCCAGGGGGGG"
probability, mismatches = base_composition_sequence_probability(pal, base_composition(seq))
assert mismatches == 0
palindrome_list.append((len(pal), seq.find(pal), probability, mismatches, pal))

#Expect this palindrome at position 0,
pal = "CCCCCCCAGGGGGGGG"
probability, mismatches = base_composition_sequence_probability(pal, base_composition(seq))
assert mismatches == 1
palindrome_list.append((len(pal), seq.find(pal), probability, mismatches, pal))

#Expect this palindrome at position 1,
pal = "CCCCCCAGGGGGGGG"
probability, mismatches = base_composition_sequence_probability(pal, base_composition(seq))
assert mismatches == 1
palindrome_list.append((len(pal), seq.find(pal), probability, mismatches, pal))

#sorted not available on python 2.3
#assert sorted(palindrome_list) == sorted(iteration_length(seq, 15,150, 1.0/len(seq)))
del seq, pal, probability, mismatches, palindrome_list

#Example Three
#=============
seq = "CCCCCCCCCCGGGGGGGGGG"
palindrome_list = iteration_length(seq, 15,150, 1.0/len(seq)) #,fh
assert palindrome_list == [(15, 2, 0.0078125, 0, 'CCCCCCCCGGGGGGG'),
                           (15, 3, 0.0078125, 0, 'CCCCCCCGGGGGGGG'),
                           (16, 1, 0.03125, 1, 'CCCCCCCCCGGGGGGG'),
                           (16, 2, 0.00390625, 0, 'CCCCCCCCGGGGGGGG'),
                           (16, 3, 0.03125, 1, 'CCCCCCCGGGGGGGGG'),
                           (17, 0, 0.03125, 1, 'CCCCCCCCCCGGGGGGG'),
                           (17, 1, 0.00390625, 0, 'CCCCCCCCCGGGGGGGG'),
                           (17, 2, 0.00390625, 0, 'CCCCCCCCGGGGGGGGG'),
                           (17, 3, 0.03125, 1, 'CCCCCCCGGGGGGGGGG'),
                           (18, 0, 0.017578125, 1, 'CCCCCCCCCCGGGGGGGG'),
                           (18, 1, 0.001953125, 0, 'CCCCCCCCCGGGGGGGGG'),
                           (18, 2, 0.017578125, 1, 'CCCCCCCCGGGGGGGGGG'),
                           (19, 0, 0.001953125, 0, 'CCCCCCCCCCGGGGGGGGG'),
                           (19, 1, 0.001953125, 0, 'CCCCCCCCCGGGGGGGGGG'),
                           (20, 0, 0.0009765625, 0, 'CCCCCCCCCCGGGGGGGGGG')]
del seq, palindrome_list
