# motif-mark

Goal: Primary goal was to be able to create an object-orientated programming script to create a visual representation of where motifs are located on all genes in a FASTA file.

Two input files:
- FASTA file: Correctly formatted FASTA file with sequences containing less than or equal to 1000 bases.
- Motif file: One motif sequence per line saved in a text file. Each motif is assumed to be less than or equal to 10 bases

Assumptions and Capabilities:
- Motifs may contain ambiguous nucleotides
- Can handle up to a maximum number of 5 motifs at a time
- Can handle overlapping motifs
- Drawing is done to scale
- Only require pycairo package to operate

Output: one well-labeled png figure per FASTA file. This figure will include all of the genes found in the original FASTA file, plus a color legend for motifs.