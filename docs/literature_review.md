# Literature Review and Theoretical Background

## 1) Dynamic programming alignment foundations

### Needleman–Wunsch (global alignment)
Needleman–Wunsch is a dynamic programming algorithm that aligns two complete sequences end-to-end.
It fills a scoring matrix from the top-left to bottom-right using recurrence relations that consider:

- match/mismatch substitution score,
- insertion (gap in sequence B),
- deletion (gap in sequence A).

The best-scoring path through the full matrix gives a global alignment. This is most appropriate when
sequences are expected to be homologous across most of their lengths.

### Smith–Waterman (local alignment)
Smith–Waterman uses a similar matrix recurrence, but resets negative scores to zero. This enables it to
find the highest-scoring local region shared by two sequences, rather than forcing an end-to-end alignment.
It is preferred when two sequences may share only a conserved domain or motif.

## 2) Substitution matrices and gap penalties

### BLOSUM/PAM matrices
Substitution matrices encode biologically plausible residue replacements.

- **BLOSUM** matrices are empirically derived from observed substitutions in conserved blocks.
  BLOSUM62 is a widely used default for protein alignments.
- **PAM** matrices model evolutionary change over accepted mutation distances.

Higher matrix scores reward likely substitutions; lower/negative scores penalize unlikely replacements.

### Gap penalties
Gap modeling typically uses an affine scheme:

- **gap-open penalty**: cost to introduce a new gap,
- **gap-extend penalty**: smaller cost to lengthen an existing gap.

This reflects the biological expectation that one longer indel is often more plausible than many short indels.

## 3) Biological significance

Sequence alignment supports key bioinformatics tasks, including:

- identifying homologous genes/proteins,
- detecting conserved catalytic residues and motifs,
- inferring functional similarity,
- supporting structure/function prediction,
- guiding evolutionary and phylogenetic analyses.

## 4) References (starter reading)

1. Needleman SB, Wunsch CD. *A general method applicable to the search for similarities in the amino acid sequence of two proteins*. J Mol Biol. 1970.
2. Smith TF, Waterman MS. *Identification of common molecular subsequences*. J Mol Biol. 1981.
3. Henikoff S, Henikoff JG. *Amino acid substitution matrices from protein blocks*. Proc Natl Acad Sci USA. 1992.
