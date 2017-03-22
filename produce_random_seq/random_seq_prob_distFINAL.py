#!/usr/bin/env python
#### Title: function to generate a random nucleotide sequence of specified
# length using a given base composition



"""
functions:      def base_composition(s)
                            output bc to be used in random_sequence_generator
                            (length, bc)
                            

                def random_sequence_generator(length, bc)
                            output seq

"""
# imports
import sys
import os
import random
from Bio.Seq import Seq
from Bio import SeqIO
from collections import Counter


if "-v" in sys.argv or "--version" in sys.argv:
    print("create random seq v0.1.0. This take a db calculates the " +
          "average ATCG and average length, then creates N number of " +
          "random sequence file")
    sys.exit(os.system(cmd))
    
############################################################################
######## function 1: This works out the codon bias of the sequence in question
########            becareful which seq it is analysing



def base_composition(s):
    """function to work out the 'C,G,A,T' content of the genome/sequecne
    in question. The function takes a string or Bio.Seq Seq as the input"""
    s1 = str(s).upper()
    #creat a empty dictionary to put the values of the bases to
    length = float (len(s1))
    bc= {'A': (float (s1.count('A')/ length)),\
         'T': float(s1.count('T')/ length),\
         'G': float(s1.count('G')/ length),\
         'C': float(s1.count('C')/ length)}
    return bc


##############################################################################
## function 2: the produces a random sequence with codon bias based on the 
##           ratio defined in bc

def random_letter_generator(length, bc={'A':0.1, 'T':0.4, 'C':0.25, 'G':0.25}):
    """Returns a single letter using the base composition as a probability
    distrubtion."""
    #create an empty seq to add the bases to
    random_seq = ""
    # while loop that stop when our new seq reaches the length defined by the
    #user. random choice picks out bases from our shuffled sequence.
    a_threshold = bc['A']
    t_threshold = a_threshold + bc['T']
    c_threshold = t_threshold + bc['C']    
    while len (random_seq) <= (length-1):
        random_number = random.random()
        # if the random generated number is less than 
        if random_number <= a_threshold:
            random_seq += 'A'
        elif random_number <= t_threshold:
            random_seq += 'T'
        elif random_number <= c_threshold:
            random_seq += 'C'
        else:
            random_seq += 'G'
    return random_seq


def write_shuffled_sequence(original_filename, shuffled_filename) :
    """Write a shuffled version of the sequence to a FASTA file."""

    record = SeqIO.read(open(original_filename), "fasta")

    data = list(record.seq)
    random.shuffle(data)
    shuffled_string = "".join(data)

    f= open(shuffled_filename, 'w')
    print >> f, ">random_version_%s" % record.id
    #Elegant approach to line wrapping, writing out line by line@
    for index in range(0, len(shuffled_string), 70) :
        print >> f, shuffled_string[index:index+70]
    f.close()

try:
    # New in Python 3.4
    from statistics import mean
except ImportError:
    def mean(list_of_values):
        """Calculate the mean average of a list of numbers."""
        # Quick and dirty, assumes already a list not an interator
        # so don't have to worry about getting the divisor.
        # Explicit float(...) to allow for Python 2 division.
        return sum(list_of_values) / float(len(list_of_values))

assert mean([1,2,3,4,5]) == 3

################# mini quick test #########################################

for length in [10, 50, 100] :
    random_seq = random_letter_generator(length)
    print random_seq, len(random_seq), "expected length", length
    print base_composition(random_seq)

############################################################################
############################## TEST CODE ###################################

#### test 1 ################################################################
    
""" this test should bring back a random sequence with length 1000 and bc as
defined below in the function"""

new_test1= random_letter_generator(1000)
print 'TEST1:\n "LENGTH 1000"',len(new_test1), base_composition(new_test1), \
"bc={'A':0.1, 'T':0.4, 'C':0.25, 'G':0.25})\n"



#### test 2###################################################################
##
##""" this test should bring back a random sequence with length 10000 and bc as
##defined below, should return no 'A' bases""" # no 'A' bases
##
##
##new_test2= random_letter_generator(10000, bc={'A':0.00, 'T':0.5, 'C':0.25,\
##                                              'G':0.25})
##print 'TEST2:\n "LENGTH 1000"',len(new_test2), base_composition(new_test2), \
##"bc={'A':0.0, 'T':0.5, 'C':0.25, 'G':0.25})\n"
##
##
####test 3####################################################################
##
##""" this test should bring back a random sequence with length 10000 and bc as
##defined below""" #should return no c or g bases
##
##new_test3= random_letter_generator(5678, bc={'A':0.5, 'T':0.5, 'C':0.0,\
##                                              'G':0.0})
##print 'TEST3:\n "LENGTH 5678"',len(new_test3), base_composition(new_test3), \
##"bc={'A':0.5, 'T':0.5, 'C':0.0, 'G':0.0})\n"
##
##
####test 4####################################################
##
##""" this test should bring back a random sequence with length 10000 and bc as
##defined below, should return even bases""" #even distribution
##
##new_test4= random_letter_generator(9999, bc={'A':0.25, 'T':0.25, 'C':0.25,\
##                                              'G':0.25})
##print 'TEST4:\n "LENGTH 9999"',len(new_test4), base_composition(new_test4), \
##"bc={'A':0.25, 'T':0.25, 'C':0.25, 'G':0.25})\n"
##
##
##
####test 5####################################################
##
##""" this test should bring back a random sequence with length 10000 and bc as
##defined below, should return g bases""" #all g bases
##
##new_test5= random_letter_generator(10000, bc={'A':0.0, 'T':0.0, 'C': 0.0,\
##                                              'G':1.0})
##print 'TEST5:\n "LENGTH 1000"',len(new_test5), base_composition(new_test5), \
##" should be 'bc={'A':0.0, 'T':0.0, 'C':0.0, 'G':1})'\n"

########
# run some code.

seq_lengths = []
# easy way of counting
base_composition_dict = Counter({'A':0.0, 'T':0.0, 'C': 0.0,'G':0.0})

final_GC = {'A':0.0, 'T':0.0, 'C': 0.0,'G':0.0}

with open ("ITS_database_NOT_confirmed_correct_last14bases_removedabundance.fasta",
           "r") as fasta_file:
    # count to decide how many to divide by later
    seq_count = 0 
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        seq_count = seq_count + 1
        # appnend value to list
        seq_lengths.append(len(seq_record.seq))
        # call function to get base composition bc
        bComp = base_composition(seq_record.seq)
        bComp = Counter(bComp)
        #print base_composition_dict
        # use the counter from collections
        base_composition_dict = base_composition_dict + bComp
        
    print "mean_lenghts : ", mean(seq_lengths)
    mean_len = mean(seq_lengths)
    print base_composition_dict

for key, vals in base_composition_dict.items():
    final_GC[key] = vals/seq_count

print final_GC

# create sequences with this codon bias using the mean len of the seq.
random_seq_num = 0
for i in range(1, 50):
    random_seq_num = random_seq_num + 1
    random_seq = random_letter_generator(mean_len, final_GC)
    #print random_seq
    name_out = "ITS_random_seq_%d.fasta" % (random_seq_num)

    f_out= open(name_out, 'w')
    out_string = ">ITS_random_vers_%d\tmean_len = %f BaseComp_av = %s\n" % (random_seq_num,
                                                                            mean_len,
                                                                            str(final_GC))
    f_out.write(out_string)
    #Elegant approach to line wrapping, writing out line by line@
    for index in range(0, len(random_seq), 70) :
        sequence = random_seq[index:index+70] + "\n"
        f_out.write(sequence)
            
    f_out.close()
    
