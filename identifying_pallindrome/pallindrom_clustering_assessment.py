#!/usr/bin/env python

""" this is a function to return start position minus the next start position
using an excisting data set.

WHY?: This will provide us with the data to see if there is any evidence of
clustering of palindromes against a random data set.I would predict the 'random'
to have a normal distribution of distance between one palindrome to the
next.... if we compare our real data against the data generated from the
random data we can prove if palindrome clustering is happening."""

#import
import os
from optparse import OptionParser
import sys
import time

#   TITLE: This is a function to return start position minus start positions
#in an ordered data set of identified palindromes. This is fo the identification
#of palindrome clustering


usage = """

Usage:

./pallindrom_clustering_assessment.py --folder default=os.getcwd() --base_name


#   TITLE: This is a function to return start position minus start positions
#in an ordered data set of identified palindromes. This is fo the identification
#of palindrome clustering
this is a function to return start position minus the next start position
using an excisting data set.

WHY?: This will provide us with the data to see if there is any evidence of
clustering of palindromes against a random data set.I would predict the 'random'
to have a normal distribution of distance between one palindrome to the
next.... if we compare our real data against the data generated from the
random data we can prove if palindrome clustering is happening.!!!
"""
parser = OptionParser(usage=usage)



parser.add_option("--folder", dest="folder", default=os.getcwd(),
                  help="previous output folder is, where the"
                  " condensed palindrome files are (real and shuffled data",
                  metavar="folder")
parser.add_option("--base_name", dest="base_name", default=None,
                  help="base_name for the species")

(options, args) = parser.parse_args()


#--folder
folder = options.folder
#--base_name
base_name = options.base_name

###########################################################################


""" function list:  load_palindromes_file(filename)
                        return palindromes (list of tuples)

                    start_length_minus_next_start_length(palindromes):
                        return
                    
                    my_outfile(outfilename, output):
                        return None

"""
###########################################################################
##function 1: this function is to open the file of choice for further work..


def load_palindromes_file(filename):
    """this is a function to open a desired file and return the contents as a
    list of tuples"""
    f= open(filename, "r")
    #assign the file contents to the variable data
    data = f.readlines()
    #remove the \n new line and \t characters
    data1 = [line.rstrip("\n").split("\t") for line in (data)
             if line.strip() != ""]
    #to convert data 1 into a list of tuples.
    #remove the title of the file
    data2 = data1 [:-1]
    #starts = [(int(s[0]),(int(s[1]))) for s in (data1)]
    #THE NEXT LINE IS SPECIFIC TO THE OVERAL TASK NOT TO THIS FUNCTION
    palindromes = [(int(s[0]), int(s[1]), float(s[2]), int(s[3]), s[4])\
                   for s in (data2)if s[0][0] != ("#")]
    #remove the title of the file
    #if there is a title
    #data2 = data1 [1:]
    f.close()
    #print data
    return palindromes


#############################################################################
# function 2: the function to do the task in hand...... return the start length
# position minus the next start length position for the further analysis with
# the goal of identifying if clustering is occuring.


def start_length_minus_next_start_length(palindromes):
    """ this is a function to return a list of numbers (integers) that are the
results of the start position minus the next start position in the data set

The data taken into this function is a list of tuples returned from the
load_palindromes_file function"""
    # to return the start positions
    starts = [(s[1]) for s in (palindromes)]
    palin = palindromes[1:]
    sorted (palin)
    #print palin
    # the sorted function orders the list for low to high
    start_length_ordered = sorted (starts)
    #print start_length_ordered
    #we want the list high to low.... so the following is used:
    longest_start_ordered = start_length_ordered[::-1]
    #print longest_start_ordered
    # empty list to append the results to
    start_length_minus_next_start_length = []
    # iteration through the start positions to return the results of the start
    # position minus the next start position in the list
    for i in range (len(longest_start_ordered)-1):
        j = i+1
        value = longest_start_ordered[i] - longest_start_ordered [j]
        start_length_minus_next_start_length.append(value)
        
    return start_length_minus_next_start_length


# function 3: function to write the results to a file

def write_start_minus_start_to_file(outfilename, output):
    """ this is a function to write the results of a program to a defined file"""
    #this opens the destined file to write 'what we want' to it ('w')
    fh= open(outfilename, "w")
    #the output needs to be formatted in a way that we want it
    #this print statement write to out outfile
    #empty string to add the formated strings to.
    output_formated2 = ""
    #for loop to convert the individual tuples(i) into formated string.
    #The formated strings need to be added to the final list...
    #then printed to the desired file.
    for i in output:
        #print 'i is : ', i
        output_formated= '%d\n' %(i)
        output_formated2+=output_formated
    #this print statement write to out outfile
    print >> fh, output_formated2
    #closes the outfile
    fh.close()
    return output_formated2 #replace with None when happy with output!

######################################################################################################################################################
## actual code......

#for loop to run the above functions on all the txt files containing the
#'condensed' palindromes.
full_list_starts = []
prefix = "%s/%s" % (folder, base_name)
random_prefix = "%s/random/%s" % (folder, base_name)
for i in range (1,31):
    start_time=time.time()
    in_filename = random_prefix+"_palin_random_condensed_N%02i.txt" % (i)
    
    data = load_palindromes_file(in_filename)

    starts = start_length_minus_next_start_length(data)
    temp_out = random_prefix+'random_distance_to_next_pallindrome_%02i.txt'%i
    write_start_minus_start_to_file(temp_out,starts)
    full_list_starts+=starts
    #print 'that took, %.3f' %(end_time - start_time)
    
result_out = random_prefix+"random_data_distance_to_next_pallindrome.txt"
write_start_minus_start_to_file(result_out,full_list_starts)
end_time=time.time()
print 'that took, %.3f' %(end_time - start_time)

################################################################################################################
#real data

in_filename = prefix+"_palin_REAL_condensed.txt"
data = load_palindromes_file(in_filename)
starts = start_length_minus_next_start_length(data)
temp_out = prefix+'REAL_distance_to_next_pallindrome.txt'
write_start_minus_start_to_file(temp_out,starts)
