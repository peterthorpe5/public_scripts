import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_nucleotide, generic_protein
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio import AlignIO

def sequence_back_translate(aligned_protein_record, unaligned_nucleotide_record):
    #TODO - Separate arguments for protein gap and nucleotide gap?
    seq = []
    nuc = str(unaligned_nucleotide_record.seq)
    for amino_acid in aligned_protein_record.seq:
        if amino_acid == "-":
            seq.append("---")
        else:
            seq.append(nuc[:3])
            nuc = nuc[3:]
    assert not nuc, "Nucleotide sequence for %r longer than protein %s" \
        % (unaligned_nucleotide_record.id, aligned_protein_record.id)

    aligned_nuc = unaligned_nucleotide_record[:] #copy for most annotation
    aligned_nuc.letter_annotation = {} #clear this
    aligned_nuc.seq = Seq("".join(seq), generic_nucleotide) #replace this
    assert len(aligned_protein_record.seq) * 3 == len(aligned_nuc)
    return aligned_nuc

def alignment_back_translate(protein_alignment, nucleotide_records):
    """Thread nucleotide sequences onto a protein alignment."""
    #TODO - Separate arguments for protein and nucleotide gap characters?
    aligned = []
    for protein in protein_alignment:
        try:
            nucleotide = nucleotide_records[protein.id]
        except KeyError:
            raise ValueError("Could not find nucleotide sequence for protein %r" \
                             % protein.id)
        aligned.append(sequence_back_translate(protein, nucleotide))
    return MultipleSeqAlignment(aligned)


align_format = "fasta" #could be clustal
prot_align_file = "p.ensembl.20150928_vs_SBP.hmm.DOMAIN_ONLY_only.aligned.refined.fasta"
nuc_fasta_file = "p.ensembl.20150928_vs_SBP.hmm_NUCLEOTIDE_domain_only.fasta"
nuc_align_file = "p.ensembl.20150928_vs_SBP.hmm.DOMAIN_ONLY_only.aligned.refined_back_trans.fasta"

prot_align = AlignIO.read(prot_align_file, align_format, alphabet=generic_protein)
nuc_dict = SeqIO.index(nuc_fasta_file, "fasta")
nuc_align = alignment_back_translate(prot_align, nuc_dict)
AlignIO.write(nuc_align, nuc_align_file, align_format)
print "Done"
