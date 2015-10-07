#!/usr/bin/env python
"""Back-translate a protein alignment to nucleotides

This tool is a short Python script (using Biopython library functions) to
load a protein alignment, and matching nucleotide FASTA file of unaligned
sequences, in order to produce a codon aware nucleotide alignment - which
can be viewed as a back translation.

The development repository for this tool is here:

* https://github.com/peterjc/pico_galaxy/tree/master/tools/align_back_trans  

This tool is available with a Galaxy wrapper from the Galaxy Tool Shed at:

* http://toolshed.g2.bx.psu.edu/view/peterjc/align_back_trans

See accompanying text file for licence details (MIT licence).

This is version 0.0.3 of the script.
"""

import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio import AlignIO
from Bio.Data.CodonTable import ambiguous_generic_by_id

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.5"
    sys.exit(0)

def stop_err(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

def check_trans(identifier, nuc, prot, table):
    """Returns nucleotide sequence if works (can remove trailing stop)"""
    if len(nuc) % 3:
        stop_err("Nucleotide sequence for %s is length %i (not a multiple of three)"
                 % (identifier, len(nuc)))

    p = str(prot).upper().replace("*", "X")
    t = str(nuc.translate(table)).upper().replace("*", "X")
    if len(t) == len(p) + 1:
        if str(nuc)[-3:].upper() in ambiguous_generic_by_id[table].stop_codons:
            #Allow this...
            t = t[:-1]
            nuc  = nuc[:-3] #edit return value
    if len(t) != len(p):
        err = ("Inconsistent lengths for %s, ungapped protein %i, "
               "tripled %i vs ungapped nucleotide %i." %
               (identifier, len(p), len(p) * 3, len(nuc)))
        if t.endswith(p):
            err += "\nThere are %i extra nucleotides at the start." % (len(t) - len(p))
        elif t.startswith(p):
            err += "\nThere are %i extra nucleotides at the end." % (len(t) - len(p))
        elif p in t:
            #TODO - Calculate and report the number to trim at start and end?
            err += "\nHowever, protein sequence found within translated nucleotides."
        elif p[1:] in t:
            err += "\nHowever, ignoring first amino acid, protein sequence found within translated nucleotides."
        stop_err(err)


    if t == p:
        return nuc
    elif p.startswith("M") and "M" + t[1:] == p:
        #Close, was there a start codon?
        if str(nuc[0:3]).upper() in ambiguous_generic_by_id[table].start_codons:
            return nuc
        else:
            stop_err("Translation check failed for %s\n"
                     "Would match if %s was a start codon (check correct table used)\n"
                     % (identifier, nuc[0:3].upper()))
    else:
        #Allow * vs X here? e.g. internal stop codons
        m = "".join("." if x==y else "!" for (x,y) in zip(p,t))
        if len(prot) < 70:
            sys.stderr.write("Protein:     %s\n" % p)
            sys.stderr.write("             %s\n" % m)
            sys.stderr.write("Translation: %s\n" % t)
        else:
            for offset in range(0, len(p), 60):
                sys.stderr.write("Protein:     %s\n" % p[offset:offset+60])
                sys.stderr.write("             %s\n" % m[offset:offset+60])
                sys.stderr.write("Translation: %s\n\n" % t[offset:offset+60])
        stop_err("Translation check failed for %s\n" % identifier)

def sequence_back_translate(aligned_protein_record, unaligned_nucleotide_record, gap, table=0):
    #TODO - Separate arguments for protein gap and nucleotide gap?
    if not gap or len(gap) != 1:
        raise ValueError("Please supply a single gap character")

    alpha = unaligned_nucleotide_record.seq.alphabet
    if hasattr(alpha, "gap_char"):
        gap_codon = alpha.gap_char * 3
        assert len(gap_codon) == 3
    else:
        from Bio.Alphabet import Gapped
        alpha = Gapped(alpha, gap)
        gap_codon = gap*3

    ungapped_protein = aligned_protein_record.seq.ungap(gap)
    ungapped_nucleotide = unaligned_nucleotide_record.seq
    if table:
        ungapped_nucleotide = check_trans(aligned_protein_record.id, ungapped_nucleotide, ungapped_protein, table)
    elif len(ungapped_protein) * 3 != len(ungapped_nucleotide):
        stop_err("Inconsistent lengths for %s, ungapped protein %i, "
                 "tripled %i vs ungapped nucleotide %i" %
                 (aligned_protein_record.id,
                  len(ungapped_protein),
                  len(ungapped_protein) * 3,
                  len(ungapped_nucleotide)))

    seq = []
    nuc = str(ungapped_nucleotide)
    for amino_acid in aligned_protein_record.seq:
        if amino_acid == gap:
            seq.append(gap_codon)
        else:
            seq.append(nuc[:3])
            nuc = nuc[3:]
    assert not nuc, "Nucleotide sequence for %r longer than protein %r" \
        % (unaligned_nucleotide_record.id, aligned_protein_record.id)

    aligned_nuc = unaligned_nucleotide_record[:] #copy for most annotation
    aligned_nuc.letter_annotation = {} #clear this
    aligned_nuc.seq = Seq("".join(seq), alpha) #replace this
    assert len(aligned_protein_record.seq) * 3 == len(aligned_nuc)
    return aligned_nuc

def alignment_back_translate(protein_alignment, nucleotide_records, key_function=None, gap=None, table=0):
    """Thread nucleotide sequences onto a protein alignment."""
    #TODO - Separate arguments for protein and nucleotide gap characters?
    if key_function is None:
        key_function = lambda x: x
    if gap is None:
        gap = "-"

    aligned = []
    for protein in protein_alignment:
        try:
            nucleotide = nucleotide_records[key_function(protein.id)]
        except KeyError:
            raise ValueError("Could not find nucleotide sequence for protein %r" \
                             % protein.id)
        aligned.append(sequence_back_translate(protein, nucleotide, gap, table))
    return MultipleSeqAlignment(aligned)

align_format = "fasta" #could be clustal
prot_align_file = "coo2_prot_aligned.fas"
nuc_fasta_file = "cds_coo2.fasta"
nuc_align_file = "coo2_prot_align_backtranslated"
table = 1

##if len(sys.argv) == 4:
##    align_format, prot_align_file, nuc_fasta_file = sys.argv[1:]
##    nuc_align_file = sys.stdout
##    table = 0
##elif len(sys.argv) == 5:
##    align_format, prot_align_file, nuc_fasta_file, nuc_align_file = sys.argv[1:]
##    table = 0
##elif len(sys.argv) == 6:
##    align_format, prot_align_file, nuc_fasta_file, nuc_align_file, table = sys.argv[1:]
##else:
##    stop_err("""This is a Python script for 'back-translating' a protein alignment,
##
##It requires three, four or five arguments:
##- alignment format (e.g. fasta, clustal),
##- aligned protein file (in specified format),
##- unaligned nucleotide file (in fasta format).
##- aligned nucleotiode output file (in same format), optional.
##- NCBI translation table (0 for none), optional
##
##The nucleotide alignment is printed to stdout if no output filename is given.
##
##Example usage:
##
##$ python align_back_trans.py fasta demo_prot_align.fasta demo_nucs.fasta demo_nuc_align.fasta
##
##Warning: If the output file already exists, it will be overwritten.
##
##This script is available with sample data and a Galaxy wrapper here:
##https://github.com/peterjc/pico_galaxy/tree/master/tools/align_back_trans
##http://toolshed.g2.bx.psu.edu/view/peterjc/align_back_trans
##""")

try:
    table = int(table)
except:
    stop_err("Bad table argument %r" % table)

prot_align = AlignIO.read(prot_align_file, align_format, alphabet=generic_protein)
nuc_dict = SeqIO.index(nuc_fasta_file, "fasta")
nuc_align = alignment_back_translate(prot_align, nuc_dict, gap="-", table=table)
AlignIO.write(nuc_align, nuc_align_file, align_format)
