# DCIT 411 — Sequence Alignment with Biopython

This repository contains a practical implementation of the **DCIT 411 Bioinformatics project**: pairwise and multiple sequence alignment using Biopython, along with benchmarking and advanced extensions.

## Project structure

- `src/sequence_alignment_project.py` — core implementation.
- `tests/test_sequence_alignment_project.py` — unit tests for core workflows.
- `examples/run_demo.py` — runnable demo using bundled sample data.
- `data/sample_sequences.fasta` — short protein sequences for experimentation.
- `docs/literature_review.md` — theoretical background and references.

## Literature review and theory (Task 1)

A concise theory section covering Needleman–Wunsch, Smith–Waterman, substitution matrices (BLOSUM/PAM),
affine gap penalties, and biological applications is provided in:

- `docs/literature_review.md`

## Data collection and preprocessing (Task 2)

Implemented in code:

- FASTA loading (`load_fasta_sequences`).
- Sequence cleaning (`preprocess_sequence`) that removes unsupported characters.

Data source note:

- `data/sample_sequences.fasta` is a small **teaching/demo dataset** for reproducible local testing.
- For formal project submission, replace/extend this with sequences downloaded from **NCBI** or **UniProt** and document accession IDs in your report.

## Pairwise alignment (Task 3)

- Global alignment (Needleman–Wunsch-style) via `PairwiseAligner(mode="global")`.
- Local alignment (Smith–Waterman-style) via `PairwiseAligner(mode="local")`.
- Matrix + gap-penalty customization.
- Quality metrics: alignment score, identity, similarity.

## Multiple sequence alignment (Task 4)

- Wrappers for ClustalW / MUSCLE / MAFFT (when installed in PATH).
- Runtime + memory benchmark helper.
- Consensus sequence generation from MSA.

### Benchmark reporting template

Populate this table after running your experiments:

| Tool     | Runtime (s) | Peak Memory (KB) | Notes on alignment quality |
|----------|-------------|------------------|----------------------------|
| ClustalW |             |                  |                            |
| MAFFT    |             |                  |                            |
| MUSCLE   |             |                  |                            |

## MSA visualization (Task 4)

This project supports Biopython alignment objects for inspection and printing.
You can additionally export alignments and inspect them in Jalview.

A simple terminal-style alignment block example:

```text
seq1  MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ
seq2  MKTAYIAKQRTISFVKSHFSRQDILDLWIFHTQGYFPQ
seq3  MKTAYIAKQRTISFVKSYFSRQDILDLWIFHTQGFFPQ
cons  MKTAYIAKQRTXSFVKSXFSRQDILDLWIFHTQGXF-PQ
```

## Advanced topics (Task 5)

- PSI-BLAST command builder utility.
- Structural alignment helper using `Bio.PDB.Superimposer`.
- Consensus sequence generation from aligned residues.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python examples/run_demo.py
```

## Testing

```bash
python -m pytest -q
```

## Notes

- External MSA tools (ClustalW/MUSCLE/MAFFT) must be installed separately for MSA execution.
- PSI-BLAST and structural alignment require local BLAST/PDB inputs in your environment.
