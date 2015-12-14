READme info for script to Identify putative HGT events
======================================================

basic usage:

``./Lateral_gene_transfer_predictor.py`` -h 



script to open up a tab blast output and generate Alien index scores,
this finds the kingdom that blast mtach has been assigned.
script: Parse the NCBI taxonomic database node.dmp to get the
taxonomic relationship

author: Peter Thorpe September 2015. The James Hutton Insitute, Dundee, UK.

./Lateral_gene_transfer_predictor.py -h
Usage: Use as follows:

$ ``./Lateral_gene_transfer_predictor.py`` -i blast_w_tax_id.tab --tax_filter_out 6656 (e.g.arthropoda) --tax_filter_up_to 33208 (e.g. metazoan) -o LTG_results.out


taxid - 6231 (nematoda)

Tax databse from NCBI is require. Download, unzip, and use -p /PATH/TO/   scripts will find them from here.

    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -zxvf taxdump.tar.gz


#6656 = filter_out_tax_id --tax_filter_out
#33208 = Metazoa   -  for me this is the tax id I want to go up to --tax_filter_up_to

What:
To determine Lateral gene transfer event (LGT). An alien index score needs to be generated. Score > 45
is a candidate LGT - however, it could be contamination. The user will have to decide.


How:

Alien index: (plagerised from Seb Eves Van Den Akker et al, Naccubus paper....)
taxonomic origin: Bacteria, Fungi, Plantae, Metazoa, and Other (including protists), and
the best expect value from each group was used to calculate an Alien Index (AI) as given
by the formula AI = log((Best E-value for Metazoa) + e-200) - log((Best E-value for Non-
Metazoa) + e-200). If no hits were found, the E-value was set to 1. Thus, AI is allowed to
vary in the interval between +460 and -460, being positive when top non-metazoan hits
yielded better E-values than the top metazoan ones. Entries with incomplete taxonomy,
such as X-ray structures, or hits belonging to the same phylum as the query sequence (
i.e Rotifera, filter_out_tax_id or Nematoda), were excluded from the analysis


    BLAST DATA should be formatted as:

qseqid = Query Seq-id (ID of your sequence)
sseqid = Subject Seq-id (ID of the database hit)
pident = Percentage of identical matches
length = Alignment length
mismatch = Number of mismatches
gapopen = Number of gap openings
qstart = Start of alignment in query
qend = End of alignment in query
sstart = Start of alignment in subject (database hit)
send = End of alignment in subject (database hit)
evalue = Expectation value (E-value)
bitscore = Bit score
salltitles = TOP description of the blast hit
staxids = tax_id
scientific_name
scomnames = common_name
sskingdoms = kingdom



MORE INFO:

Alien index:  (http://www.sciencemag.org/content/suppl/2008/05/29/320.5880.1210.DC1/Gladyshev.SOM.pdf)
    Massive Horizontal Gene Transfer in Bdelloid Rotifers :
    taxonomic origin: Bacteria, Fungi, Plantae, Metazoa, and Other (including protists), and
    the best expect value from each group was used to calculate an Alien Index (AI) as given
    by the formula AI = log((Best E-value for Metazoa) + e-200) - log((Best E-value for Non-
    Metazoa) + e-200). If no hits were found, the E-value was set to 1. Thus, AI is allowed to
    vary in the interval between +460 and -460,being positive when top non-metazoan hits
    yielded better E-values than the top metazoan ones. Entries with incomplete taxonomy,
    such as X-ray structures, or hits belonging to the same phylum as the query sequence (
    i.e Rotifera, filter_out_tax_id or Nematoda), were excluded from the analysis

    Based on AI, genes were classified as foreign (AI>45), indeterminate (0<AI<45), or metazoan (AI<0)

There will be Four outfiles, using the -o prefix.

1) prefix name :
 This will contain the best metazoan and non metazoan hit per blast query. Depending on your filters specified.

2) prefix_precursor_value_temp :
 this file contain the precursor value for each of those if out_file(1). This is used for the next file

3) prefix_Alien_index.out :
This file contains all the AI scores for each BLAST query. In one long line the best metazoanan
and non-metazoan hits can be seen. A final comment is added if the programs believes there to be a HGT/LGT event,
or if this is contamination.

Contamination is based on a high AI score and greater than 70% identity to a non-metazoan. This can ofcourse be changed to suit the user.

4) Perfix_LGT_candifates.out :
This file contains all scores greater than 0. The final comments box is a note to say if it think it is potential
contamination or if it may be a HGT/LTG event.



Options:
  -h, --help            show this help message and exit
  -i FILE, --in=FILE    the tab output from blast/diamond. This must have
                        tax_id info in this!!If you use diamond, please get
                        this info using 'add_taxonomic_info_to_tab_output.py'
  -p PATH, --path=PATH  Directory containing relevant taxonomy/database files
                        Default is the current working directory. This is not
                        used with the main input and output filenames.
  --tax_filter_out=TAX_FILTER_OUT
                        The tax ID to filter out: for this analysis the Phylum
                        which your BEASTof interest if found. e.g. Aphids are
                        from Arthropoda, therefore this would be 6656, whihc
                        is the dwefault value. This will filter out all blast
                        hit which are from this phylum. It is possible to put
                        a species/kingdom tax_id in here ... whatever floats
                        your boat.
  --tax_filter_up_to=TAX_FILTER_UP_TO
                         The tax_id to 'walk up to', to determine assignment.
                        By default this is metazoa.The script work out the
                        best metazoan to non-metazoan hit. But this can be
                        altered if you wish to alter this
  --tax_coloumn=TAX_COLOUMN
                        the coloumn with the tax_id info. Defulat is 14(as
                        counted by a human/ not a computer
  -o FILE, --out=FILE   Output filename - default=
                        infile__tab_blast_LGT_results
