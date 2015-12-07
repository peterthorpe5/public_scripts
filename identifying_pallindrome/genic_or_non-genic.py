#!/usr/bin/env python
""" this is a function to see if the 'RANROM identified palindromes'
are intergenic or
genic"""

#import
import os
from optparse import OptionParser
import sys
import time

#############################################################################
#MAIN TITLE: PALINDROME INTERGENIC OR NOT?


#   TITLE: This is a funtion to compare the results of the indentified palindromes
#against the starts and end positions of the genome in question.

#why?: if the palindrome occurs within the boundries of a gene start....end
#position, then this is said to be a genic palindrome. Palindromes are expected
#to be found in intergenic regions, not within the boundires of the gene.
#The result from this code will help us with the hypothesis, with the comparison
#using the random data as well....!!!

###########################################################################

usage = """

Usage:

./genic_or_non-genic.py --ptt  --gff --folder (default=os.getcwd()) --base_name

#MAIN TITLE: PALINDROME INTERGENIC OR NOT?


#   TITLE: This is a funtion to compare the results of the indentified palindromes
#against the starts and end positions of the genome in question.

#why?: if the palindrome occurs within the boundries of a gene start....end
#position, then this is said to be a genic palindrome. Palindromes are expected
#to be found in intergenic regions, not within the boundires of the gene.
#The result from this code will help us with the hypothesis, with the comparison
#using the random data as well....!!!
"""
parser = OptionParser(usage=usage)

parser.add_option("--ptt", dest="ptt", default=None,
                  help="coordintes for gene - can leave this option blank",
                  metavar="file")
parser.add_option("--gff", dest="gff", default=None,
                  help="gff3 file for gene predictions.",
                  metavar="file")
parser.add_option("--folder", dest="folder", default=os.getcwd(),
                  help="previous output folder is, where the"
                  " condensed palindrome files are (real and shuffled data",
                  metavar="folder")
parser.add_option("--base_name", dest="base_name", default=None,
                  help="base_name for the species")

(options, args) = parser.parse_args()

#--ptt
ptt = options.ptt
#--gff
gff = options.gff
#--folder
folder = options.folder
#--base_name
base_name = options.base_name

########################################################################################
""" function list:  load_palindromes_file(filename):
                        return palindromes (list of tuples)


                    load_gene_starts_stops_pttfile(filename):
                        return gene_start_stops
                   
                        

                    palindrome_intergenic_tester(palindromes, starts_stops):
                        return real_palindrome_list

                    
                    my_outfile(outfilename, output):
                        return None

"""

"""
structure of script:

1)load in the indentified palindromes from file
2) load in the start and stops for the gene in the genome in question.
3)Do the palindromes 'fit' inside any of the genes. The word fit can be
defined as:

    start of palindrome > start of gene AND
    end of palindrome < end of gene.
    
The identified palindrome that fit inside or do not fit inside genes need
to be scored and reurned in their respective lists.
"""

###########################################################################
##function 1: this function is to open the file of choice for further work..
# in this script the palindromes


def load_palindromes_file(filename):
    """this is a function to open a desired file and return the contents as a
    list of tuples"""
    #####!!!!!!! CHECK THIS WORKS, NEW 27TH MARCH!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #assign the file contents to the variable data
##    data = [line.strip for line in open(filename, 'rU').readlines if
##    line.strip() != ""]
    f= open(filename, "ru")
    #assign the file contents to the variable data
    data = f.readlines()
    #remove the \n new line and \t characters
    data1b = [line.rstrip("\n").rstrip().split("\t") for line in data
              if line.strip() != ""]
    #to convert data 1 into a list of tuples.
    #remove the title of the file
    #data2 = data1a [1:]
    #THE NEXT LINE IS SPECIFIC TO THE OVERAL TASK NOT TO THIS FUNCTION
    palindromes = [(int(s[0]), int(s[1]), float(s[2]), int(s[3]), s[4])\
                   for s in data1b if s[0][0] != ("#")]  #if s!="\n"
    #remove the title of the file
    #if there is a title
    #data2 = data1 [1:]
    f.close()
    #print palindromes
    return palindromes

##############################################################################
# function 2: function to load the .ptt file which contains the start and end
#positions for each gene.

def load_gene_starts_stops_pttfile(filename):
    """this is a function to open a desired file and return the contents as a
    .list of tuples containing two integers"""
    fh= open(filename, "r")
    data = fh.readlines()
    #remove the \n new line and \t characters
    #the first 3 lines are headers, so need to be removed, hence [3:]
    list_of_strings = [line.rstrip("\n").split("\t")[0] for line in data[3:]]
    #remove the .. between the values of interest...
    list_of_strings_info_split = [line.split("..") for line in (list_of_strings)]
    #the values we are interested in start and stop respectivly are in the form
    #of strings currently. we need to convert these to tuple. Therefore, a list
    #containing tuple.
    #to convert list_of_strings_info_split into a list of tuples.
    gene_start_stops = [(int(i[0]), int(i[1])) for i in\
                        list_of_strings_info_split]
    #equivelant for lopp for the above list comprehension
##    for i in list_of_strings_info_split:
##        start_stops = int(i[0]), int(i[1])
##        gene_start_stops.append(start_stops)
    #print 'GENE: STARTS, STOPS', gene_start_stops
    return gene_start_stops

######################################################################################
# function to load from GFF3
def load_gene_starts_stops_GFFfile(filename):
    """this is a function to open a desired file and return the contents as a
    .list of tuples containing two integers"""
    fh= open(filename, "r")
    data = fh.readlines()
    #remove the \n new line and \t characters
    #the first 3 lines are headers, so need to be removed, hence [3:]
    list_of_strings = [line.rstrip("\n").split("\t") for line in data if not line.startswith("#")]
    #remove the .. between the values of interest...
    #list_of_strings_info_split = [line.split("..") for line in (list_of_strings)]
    #the values we are interested in start and stop respectivly are in the form
    #of strings currently. we need to convert these to tuple. Therefore, a list
    #containing tuple.
    #to convert list_of_strings_info_split into a list of tuples.
    #print list_of_strings
    gene_start_stops = [(int(i[3]), int(i[4])) for i in\
                        list_of_strings]

    return gene_start_stops


#############################################################################
# function 3: the function to do the task in hand...... remove repeated palindromes

def palindrome_intergenic_tester(palindromes, starts_stops):
    """this is a function to compare the starts and end of the plaindromes
    from the imported file and compare these to the starts and stops for the
    genes for the genome in question

    Do the palindromes 'fit' inside any of the genes. The word fit can be
    defined as:

    start of palindrome > start of gene AND
    end of palindrome < end of gene.
    
    The identified palindrome that fit inside or do not fit inside genes need
    to be scored and reurned in their respective lists.
    """
    
    """
    a palindrome should be: [(length, start, probability, mismatches, 'seq')]
    starts_stops should be [(int, int), (int,int)()() etc] list of tuples
    containing two integers.
    """
    #need to iterate through the palindromes list and compare to the starts
    #and stops in the genome data list.
    intragenic_score = 0
    # intragenic = within a gene
    intragenic_list = []
    intergenic_score = 0
    intergenic_list = []
    for x in palindromes:
        #print 'start of main loop \n'
        #print 'Putative palindrome x: ', x
        start= x[1]
        length= x[0]
        end=start+length
        #print start, end
        #we need to compare the potential palindrome to the start and end of the
        #genes in the ptt file. To do this we need to iterate
        #through the palindromes list and compare the start and end of the
        #genes_start_stops list.
        for i in starts_stops:
            #print 'Gene i: ', i
            #print i
            start_gene = i[0]
            end_gene = i[1] 
            #this if statement should test to see if the test palindrome fits
            #inside any of our 'Genes'.
            #if it does, add this result to intragenic... break the loop.
            #print 'starts: is', start, '>', 'gene start', start_gene
            #print 'ends: is', end, '<', 'gene end', end_gene
            if start >= start_gene and end <= end_gene:
                if x not in intragenic_list:
                    intragenic_list += [x]
                    intragenic_score += 1                
                    intergenic = False
                    #print 'x is in a Gene \n'
                    #print 'I AM PUTTING X IN INTRAGENIC LIST '
                    #print 'intragenic = within a gene\n'
                break
            else:
                #this if statement should test to tsee if the potential
                #palindrome is already contained in our 'real/unique list'
                #this is a last test to reduce any potential duplicate sequences
                #print 'got to here...'
                intergenic = True
        if intergenic:
            if x not in intergenic_list:
                #print len(real_palindrome_list), 'length aplindrome list'
                pal2 = [x]
                #print 'pal2', pal2
                #print 'x is not in a gene'
                #print 'I AM PUTTING x IN intergenic list...!!! \n'
                intergenic_list += pal2
                intergenic_score += 1
                #print 'end of loop \n'      
    return intergenic_list, intergenic_score, intragenic_list, intragenic_score
    

#############################################################################
# function 3: function to write the results to a file

def write_results_to_file(outfilename, output):
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
    #return output_formated2 #replace with None when happy with output!


def write_real_results_to_file(outfilename1, outfilename2, outfilename3,
                          outfilename4, data):
    """ this is a function to write the results of a program to a defined file"""
    #this opens the destined file to write 'what we want' to it ('w')
    fh1= open(outfilename1, "w")
    fh2 = open(outfilename2, "w")
    fh3 = open(outfilename3, "w")
    fh4 = open(outfilename4, "w")
    #the output needs to be formatted in a way that we want it
    #this print statement write to out outfile
    #empty string to add the formated strings to.
    intergenic_list_formated2 = ""
    intragenic_list_formated2 = ""
    #for loop to convert the individual tuples(i) into formated string.
    #The formated strings need to be added to the final list...
    #then printed to the desired file.
    intergenic_list = data[0]
    intergenic_score = str (data [1])
    intragenic_list = data[2]
    intragenic_score = str (data[3])
    for i in intergenic_list:
        output_formated= '%d\t%d\t%.2e\t%d\t%s\n' %(i[0], i[1],i[2],i[3],i[4])
        intergenic_list_formated2+=output_formated
    for i in intragenic_list:
        output_formated= '%d\t%d\t%.2e\t%d\t%s\n' %(i[0], i[1],i[2],i[3],i[4])
        intragenic_list_formated2+=output_formated
    #this print statement write to out outfile

    fh1.write(intergenic_list_formated2)
    fh2.write(intergenic_score)
    fh3.write(intragenic_list_formated2)
    fh4.write(intragenic_score)    
    #Using print >> file will add a newline at the end!
    #print >> fh, output_formated2

    #closes the outfile
    fh1.close()
    fh2.close()
    fh2.close()
    fh4.close()
    return None

#############################################################################



#############################################################################
## actual code......

#for loop to run the above functions on all the txt files containing the
#'uncondensed' palindromes.

start_time=time.time()
intergenic_score=0
intergenic_score_list = []
intragenic_score=0
intragenic_score_list = []

if not ptt == None:
    genes_start_stops = load_gene_starts_stops_pttfile(ptt)


if not gff == None:
    genes_start_stops = load_gene_starts_stops_GFFfile(gff)
else:
    print ("no files - GFF or ptt!!")


######################################################################
#run it for the real data
filname = "%s/%s_palin_REAL_condensed.txt" \
                     % (folder, base_name)
palindromes = load_palindromes_file(filname)
        

results = palindrome_intergenic_tester(palindromes,genes_start_stops)

#print 'should be adding', results[1], 'to intergenic_score'
intergenic_score+= results[1]
intergenic_score_list.append(results[1])
#print '..and', results[3], 'to intragenic_score'
intragenic_score+=results[3]
intragenic_score_list.append(results[3])
           
out_file_prefix = "%s/%s" \
                     % (folder, base_name)

write_real_results_to_file(out_file_prefix+'real_INTERgenic_list.txt',\
                            out_file_prefix+'real_intergenic_score.txt',\
                            out_file_prefix+'real_INTRAgenic_list.txt', \
                            out_file_prefix+'real_intragenic_score.txt',\
                          results)
#write_results_to_file(out_file_prefix+'REAL_intRAgenic_full_list.txt', intragenic_score_list)
#write_results_to_file(out_file_prefix+'REAL_intERgenic_full_list.txt', intergenic_score_list)



######################################################################
#run it for the shuffled data
intergenic_score=0
intergenic_score_list = []
intragenic_score=0
intragenic_score_list = []
for i in range (1,31):
    filname = "%s/random/%s_palin_random_condensed_N%02i.txt" \
                     % (folder, base_name,i)
    palindromes = load_palindromes_file(filname)
    results = palindrome_intergenic_tester(palindromes,genes_start_stops)

    #print 'i = ', i
    #print 'should be adding', results[1], 'to intergenic_score'
    intergenic_score+= results[1]
    intergenic_score_list.append(results[1])
    #print '..and', results[3], 'to intragenic_score'
    intragenic_score+=results[3]
    intragenic_score_list.append(results[3])

out_file_prefix = "%s/random/%s" \
                     % (folder, base_name)

write_results_to_file(out_file_prefix+'random_intRAgenic_full_list.txt', intragenic_score_list)
write_results_to_file(out_file_prefix+'random_intERgenic_full_list.txt', intergenic_score_list)
end_time=time.time()
print 'that took, %.3f' %(end_time - start_time)

