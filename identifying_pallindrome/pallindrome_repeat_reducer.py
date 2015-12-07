#!/usr/bin/env python

""" this is a function to reduce the number of repeated palindromes found"""

#import
import os
from optparse import OptionParser
import sys
import time

#   TITLE: this is a function to reduce the number of repeated palindromes found
#           in a txt file. Repeated palindromes are: Palindromes that 'fit'
#           inside another palindrome sequence. For example p and m. m is the longer
#           sequence, if p start > m start and p end <  m end, p fits in m. From our
#           results file we need to reduce palindromes that 'fit' inside another...



usage = """
Usage:

./pallindrome_repeat_reducer.py --folder /path/to/ (default= os.getcwd()) --base_name (name_for_species)

#   TITLE: this is a function to reduce the number of repeated palindromes found
#           in a txt file. Repeated palindromes are: Palindromes that 'fit'
#           inside another palindrome sequence. For example p and m. m is the longer
#           sequence, if p start > m start and p end <  m end, p fits in m. From our
#           results file we need to reduce palindromes that 'fit' inside another...

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


""" function list:  my_infile(filename)
                        return palindromes (list of tuples)
                        

                    repeated_palindrome_searcher(filename):
                        return real_palindrome_list

                    
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

    #THE NEXT LINE IS SPECIFIC TO THE OVERAL TASK NOT TO THIS FUNCTION
    palindromes = [(int(s[0]), int(s[1]), float(s[2]), int(s[3]), s[4])\
                   for s in (data1) if s[0][0] != ("#")]
    #remove the title of the file
    #if there is a title
    #data2 = data1 [1:]
    f.close()
    #print data
    return palindromes 


#############################################################################
# function 2: the function to do the task in hand...... remove repeated palindromes


def repeated_palindrome_searcher(palindromes):
    """this is a function to search through a data set (list of tuple - returned
    from my_infile(filename) function). If for example, p and m. M is the longer
    sequence, if p start > m start and p end <  m end, p fits in m. From our
    result file we need to reduce palindromes that 'fit' inside another. The
    output file should be a new data set containing only unique palindrome derived
    from the content of the input file"""
    
    """ a palindrome should be: [(length, start, probability, mismatches, 'seq')]
    """
    #if the list is not ordered by length then this will ensure that they are for
    #purpose of the rest of the function. Having the list ordered is length longest
    #to smallest is ensential for the logic of the script...
    ordered_palindrome = sorted(palindromes)
    #[::-1] actually reverses the list, now they are in order, longest to smallest
    #length
    longest_first_original_palindromes = ordered_palindrome [::-1]
    #a new list to put unique plaindromes into
    real_palindrome_list = []
    #the longest palindrome cannot fit in any other sequence to this is put into the
    #list
    real_palindrome_list.append(longest_first_original_palindromes[0])
    #for loop to iterate over the tuples in the list of longest_first_original_palindromes
    #we are interested in the start and end position of each potential palidrome,
    #this will allow us to see if the start and end fits within another identified
    #palindrome.    
    for x in longest_first_original_palindromes:
        start= x[1]
        length= x[0]
        end=start+length
        #print start, end
        #we need to compare the potential palindrome to the start and end of the
        #identified 'real/unique' palindromes. To do this we need to iterate
        #through the real_palindrome_list and compare the start and end of the
        #potential and thereal palindromes to see if the potential
        #palindrome is unique.
        for i in real_palindrome_list:
            #print 'i: ', i
            #print i
            start_unique = i[1]
            length_unique = i[0]
            end_unique = start_unique+length_unique
            #this if statement should test to see if the test palindrome fits
            #inside any of our identified 'real/unique' palindromes.
            #if it does, the break the loop, if not add to real/unique list
            if start >= start_unique and end <= end_unique:
                unique_palindrome = False
                #print 'x NOT UNIQUE \n'
                break
            else:
                #this if statement should test to tsee if the potential
                #palindrome is already contained in our 'real/unique list'
                #this is a last test to reduce any potential duplicate sequences
                #print 'got to here...'
                unique_palindrome = True
        if unique_palindrome:
            if x not in real_palindrome_list:
                #print len(real_palindrome_list), 'length aplindrome list'
                pal2 = [x]
                #print 'pal2', pal2
                #print 'x is UNIQUE'
                #print 'I AM PUTTING x IN real_palindrome_list...!!!'
                real_palindrome_list+=pal2
                #print 'end of loop \n'
    return real_palindrome_list
    
#############################################################################
# function 3: function to write the results to a file

def write_palindromes_to_file(outfilename, output):
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
    print >> fh,"#length\tposition\tprobability\tmismatches\tpalindrome_seq"
    for i in output:
        output_formated= '%d\t%d\t%.2e\t%d\t%s\n' %(i[0], i[1],i[2],i[3],i[4])
        output_formated2+=output_formated
    #this print statement write to out outfile
    print >> fh, output_formated2
    #closes the outfile
    fh.close()
    return output_formated2 #replace with None when happy with output!

#############################################################################

#############################################################################
## actual code......

#for loop to run the above functions on all the txt files containing the
#'uncondensed' palindromes.
#full_true_palindrome_list = []



#for the real data
file_to_open ="%s/%s_palin_real.txt" \
                     % (folder, base_name)

data = load_palindromes_file(file_to_open)
true_palindrome_list = repeated_palindrome_searcher(data)
#full_true_palindrome_list+=true_palindrome_list
file_out = '%s/%s_palin_REAL_condensed.txt' %(folder, base_name)
write_palindromes_to_file(file_out,true_palindrome_list)
#print true_palindrome_list
print len(true_palindrome_list)


# for the shuffled data
for i in range (1,31):
    start_time=time.time()
    file_to_open = "%s/random/%s_palin_random_N%02i.txt" \
                     % (folder, base_name,i)
    data = load_palindromes_file(file_to_open)
    true_palindrome_list = repeated_palindrome_searcher(data)
    #full_true_palindrome_list+=true_palindrome_list
    out_file_name = '%s/random/%s_palin_random_condensed_N%02i.txt' %(folder, base_name,i)
    write_palindromes_to_file(out_file_name,true_palindrome_list)
#print true_palindrome_list
print len(true_palindrome_list)
end_time=time.time()
print 'that took, %.3f' %(end_time - start_time)

