
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def name_set(in_names):
    name_set = set([])
    for i in in_names:
        if not i.startswith("#"):
            name_set.add(i)
    return name_set
    

def reformat_as_fasta(filename,lenth,wanted,not_wanted, outfile):
    "this function re-write a file as a fasta file"
    f= open(outfile, 'w')
    if wanted:
        names = wanted.readlines()
        wanted_data = [line.replace("\t", "").rstrip("\n") for line in names
                  if line.strip() != ""]
        wanted_name_set = name_set(wanted_data)

    if not_wanted:
        names = not_wanted.readlines()
        not_wanted_data = [line.replace("\t", "").rstrip("\n") for line in names
                  if line.strip() != ""]
        not_wanted_name_set = name_set(not_wanted_data)
        
    #print wanted_data
    for seq_record in SeqIO.parse(filename, "fasta"):
        if len(seq_record.seq)> int(lenth):
            if wanted:
                if seq_record.id in wanted_name_set:
                    SeqIO.write(seq_record, f, "fasta")
            if not_wanted:
                if not seq_record.id in not_wanted_name_set:
                    SeqIO.write(seq_record, f, "fasta")
            else:
                SeqIO.write(seq_record, f, "fasta")                    
    
    f.close()
    return True



if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

converts

$ python rewrite_as_fasta.py -i in.fasta -l min_length_of_seq (default(3)) --not_wanted --wanted -o out.fasta

script either reformats badly formated fasta file. Within reason. Ie. Word documents will still break it.

if lenght of seq is longer than -l default 3 - writes to file.

can filter fasta by giving it a list of --not_wanted or --wanted names.

"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="in_file", default=None,
                  help="current fasta you want to reformat")

parser.add_option("--wanted", dest="wanted", default=False,
                  help="a text file with a list of wanted gene names found in the fasta file",
                  metavar="FILE")

parser.add_option("--not_wanted", dest="not_wanted", default=False,
                  help="a text file with a list of not_wanted gene names found in the fasta file",
                  metavar="FILE")

parser.add_option("-l", "--lenth", dest="lenth", default="3",
                  help="Output filename",
                  metavar="FILE")
parser.add_option("-o", "--out", dest="out", default=None,
                  help="Output filename",
                  metavar="FILE")




(options, args) = parser.parse_args()

in_file = options.in_file
not_wanted = options.not_wanted
wanted = options.wanted
out = options.out
lenth = options.lenth
convert_files(in_file, out)


reformat_as_fasta(in_file,lenth,wanted,not_wanted, out)
print 'done'

