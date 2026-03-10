# DCIT 411 — Sequence Alignment with Biopython

This repository contains a practical implementation of the **DCIT 411 Bioinformatics project**: pairwise and multiple sequence alignment using Biopython, along with benchmarking and advanced extensions.

## Project structure

- `src/sequence_alignment_project.py` — core implementation.
- `tests/test_sequence_alignment_project.py` — unit tests for core workflows.
- `examples/run_demo.py` — runnable demo using bundled sample data.
- `data/sample_sequences.fasta` — short protein sequences for experimentation.

## Features implemented

1. **Data collection and preprocessing**
   - FASTA loading.
   - Sequence cleaning (remove gaps/unknowns, uppercase conversion).

2. **Pairwise alignment**
   - Global alignment (Needleman–Wunsch-style) via `PairwiseAligner(mode="global")`.
   - Local alignment (Smith–Waterman-style) via `PairwiseAligner(mode="local")`.
   - Matrix + gap penalty customization.
   - Metrics: score, identity, similarity.

3. **Multiple sequence alignment (MSA)**
   - Wrapper for ClustalW / MUSCLE / MAFFT (when installed in PATH).
   - Runtime + memory benchmark helper.
   - MSA consensus sequence generation.

4. **Advanced topics support**
   - PSI-BLAST command builder utility.
   - Basic structural alignment helper using `Bio.PDB.Superimposer`.

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
- PSI-BLAST/structural alignment helpers are provided as practical hooks for project extensions.
