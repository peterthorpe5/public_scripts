#Title: Python script to draw the genome using, genome diagram

"""
this is an untaural process of voodoo to make a very pretty picture!!!

"""

############################################################################
#imports.......

from reportlab.lib import colors
from reportlab.lib.units import cm
# Biopython core
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Bio.Graphics.GenomeDiagram
from Bio.Graphics.GenomeDiagram import Diagram

##############################################################################
### function 1: load the palindromes start and work out the end position


def load_palindromes_file(filename):
    """this is a function to open a desired file and return the contents as a
    list of tuples"""
    f= open(filename, "r")
    #assign the file contents to the variable data
    data = f.readlines()
    #remove the \n new line and \t characters
    data1 = [line.rstrip("\n").split("\t") for line in (data)]
    #to convert data 1 into a list of tuples.
    #remove the title of the file
    data2 = data1 [:]
    print data1 [-1]
    #starts = [(int(s[0]),(int(s[1]))) for s in (data1)]
    #THE NEXT LINE IS SPECIFIC TO THE OVERAL TASK NOT TO THIS FUNCTION
    palindromes = [(int(s[0]), int(s[1]), float(s[2]), int(s[3]), s[4])\
                   for s in (data2)]
    #remove the title of the file
    #if there is a title
    #data2 = data1 [1:]
    f.close()
    #print data
    #we are only interested in the start and end of the palindromes...
    palindromes_list = []
    for i in palindromes:
        start= i[1]
        length= i[0]
        end=start+length
        info = (start,end)
        palindromes_list.append(info)
        #palindrome_list should be a list of tuples containing the start and
        #end for each palindrome
    return palindromes_list




##############################################################################
#call the function to get the 'palindromes' assigned toa variable to futher use

#intragenic - currently set to the conserved region
palindromes_intragenic = load_palindromes_file('Gr1real_INTRAgenic_list.txt')
palindromes_intergenic = load_palindromes_file('Gr1real_INTERgenic_list.txt')
#print palindromes

################################################################################

#load the genbank file that contains the genes
gbk_filename = "Gr1.gbk"
genbank_entry = SeqIO.read(open(gbk_filename), "genbank")

gdd = Diagram('Test Diagram')

#Add a track of features,
gdt_features = gdd.new_track(1, greytrack=True,
                             name="CDS Features",
                             scale_largetick_interval=10000,
                             scale_smalltick_interval=1000,
                             scale_fontsize=4,
                             scale_format = "SInt",
                             greytrack_labels=True, #e.g. 5
                             height=0.75)

#We'll just use one feature set for these features,
gds_features = gdt_features.new_set()


#And add some palindromes...
for start, end in palindromes_intragenic :
    feature = SeqFeature(FeatureLocation(start, end)) #Ideally in python 0 based
    #leaving sigil="BOX" as default
    #leaving strand as both by default
    gds_features.add_feature(feature, color="black",
                             label=False)

#And add some palindromes...
for start, end in palindromes_intergenic :
    feature = SeqFeature(FeatureLocation(start, end)) #Ideally in python 0 based
    #leaving sigil="BOX" as default
    #leaving strand as both by default
    gds_features.add_feature(feature, color="yellow",
                             label=False)
    
for feature in genbank_entry.features:
    #if feature.type not in ["CDS", "tRNA", "rRNA"] :
    if feature.type in ["source", "gene"] :
        #We're going to ignore these (ignore genes as the CDS is enough)
        continue

    #Note that I am using strings for color names, instead
    #of passing in color objects.  This should also work!
    if feature.type == "tRNA" :
        color = "red"
    elif feature.type == "rRNA":
        color = "purple"
    elif feature.type != "CDS" :
        color = "lightgreen"
    elif len(gds_features) % 2 == 0 :
        color = "blue"
    else :
        color = "green"

    gds_features.add_feature(feature, color=color,
                             sigil="ARROW",
                             arrowshaft_height=0.6,
                             arrowhead_length=0.5,
                             label_position = "start",
                             label_size = 4,
                             label_angle = 90,
                             label=True)


    

#And draw it...
gdd.draw(format='linear', orientation='landscape',
         tracklines=False, pagesize='A3', fragments=10)
gdd.write("Gr1.gbk_linear_codingPALINbox.pdf", 'PDF')
#gdd.write("NC_005213_linear.svg", 'SVG')

#And a circular version
#Change the order and leave an empty space in the center:
gdd.move_track(1,3)
gdd.draw(format='circular', tracklines=False, pagesize=(30*cm,30*cm))
gdd.write("Gr1.gbk_circular_codingPALINbox.pdf", 'PDF')
#gdd.write("NC_005213_circular.svg", 'SVG')

print "Done"
