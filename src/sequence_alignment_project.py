"""Sequence alignment utilities for DCIT 411 project."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import tracemalloc
from typing import TYPE_CHECKING, Iterable, Sequence

if TYPE_CHECKING:
    from Bio import Align
    from Bio.Align import Alignment, MultipleSeqAlignment
    from Bio.SeqRecord import SeqRecord


@dataclass
class PairwiseAlignmentResult:
    score: float
    aligned_seq_a: str
    aligned_seq_b: str
    identity: float
    similarity: float


@dataclass
class BenchmarkResult:
    method: str
    elapsed_seconds: float
    peak_memory_kb: float
    n_sequences: int


def preprocess_sequence(raw_sequence: str) -> str:
    """Normalize protein sequence by removing unsupported characters and uppercasing."""
    allowed = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
    cleaned = [ch for ch in raw_sequence.upper() if ch in allowed]
    return "".join(cleaned)


def load_fasta_sequences(path: str | Path) -> list["SeqRecord"]:
    """Load FASTA records and preprocess sequence symbols."""
    fasta_path = Path(path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file does not exist: {fasta_path}")

    from Bio import SeqIO
    from Bio.Seq import Seq

    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in: {fasta_path}")

    for record in records:
        record.seq = Seq(preprocess_sequence(str(record.seq)))
    return records


def _build_aligner(
    mode: str,
    matrix_name: str = "BLOSUM62",
    gap_open: float = -10,
    gap_extend: float = -0.5,
) -> "Align.PairwiseAligner":
    from Bio import Align
    from Bio.Align import substitution_matrices

    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = substitution_matrices.load(matrix_name)
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    return aligner


def _gapped_strings(alignment: "Alignment") -> tuple[str, str]:
    """Build gapped alignment strings from an Alignment object."""
    seq_a, seq_b = str(alignment.sequences[0]), str(alignment.sequences[1])
    coords = alignment.coordinates

    gapped_a: list[str] = []
    gapped_b: list[str] = []

    for i in range(coords.shape[1] - 1):
        a0, a1 = int(coords[0, i]), int(coords[0, i + 1])
        b0, b1 = int(coords[1, i]), int(coords[1, i + 1])

        a_chunk = seq_a[a0:a1]
        b_chunk = seq_b[b0:b1]

        len_a = a1 - a0
        len_b = b1 - b0

        if len_a == len_b:
            gapped_a.append(a_chunk)
            gapped_b.append(b_chunk)
        elif len_a > 0 and len_b == 0:
            gapped_a.append(a_chunk)
            gapped_b.append("-" * len_a)
        elif len_b > 0 and len_a == 0:
            gapped_a.append("-" * len_b)
            gapped_b.append(b_chunk)

    return "".join(gapped_a), "".join(gapped_b)


def _identity_similarity(
    aligned_a: str,
    aligned_b: str,
    matrix_name: str,
) -> tuple[float, float]:
    from Bio.Align import substitution_matrices

    matrix = substitution_matrices.load(matrix_name)
    matches = 0
    similars = 0
    positions = 0

    for a, b in zip(aligned_a, aligned_b):
        if a == "-" or b == "-":
            continue
        positions += 1
        if a == b:
            matches += 1
        if matrix.get((a, b), matrix.get((b, a), -1)) > 0:
            similars += 1

    if positions == 0:
        return 0.0, 0.0

    return matches / positions, similars / positions


def _pairwise_alignment(
    mode: str,
    seq_a: str,
    seq_b: str,
    matrix_name: str = "BLOSUM62",
    gap_open: float = -10,
    gap_extend: float = -0.5,
) -> PairwiseAlignmentResult:
    seq_a_clean = preprocess_sequence(seq_a)
    seq_b_clean = preprocess_sequence(seq_b)
    if not seq_a_clean or not seq_b_clean:
        raise ValueError("Sequences must contain at least one valid residue after preprocessing")

    aligner = _build_aligner(mode, matrix_name, gap_open, gap_extend)
    alignment = aligner.align(seq_a_clean, seq_b_clean)[0]
    aligned_a, aligned_b = _gapped_strings(alignment)
    identity, similarity = _identity_similarity(aligned_a, aligned_b, matrix_name)

    return PairwiseAlignmentResult(
        score=float(alignment.score),
        aligned_seq_a=aligned_a,
        aligned_seq_b=aligned_b,
        identity=identity,
        similarity=similarity,
    )


def pairwise_global_alignment(
    seq_a: str,
    seq_b: str,
    matrix_name: str = "BLOSUM62",
    gap_open: float = -10,
    gap_extend: float = -0.5,
) -> PairwiseAlignmentResult:
    return _pairwise_alignment("global", seq_a, seq_b, matrix_name, gap_open, gap_extend)


def pairwise_local_alignment(
    seq_a: str,
    seq_b: str,
    matrix_name: str = "BLOSUM62",
    gap_open: float = -10,
    gap_extend: float = -0.5,
) -> PairwiseAlignmentResult:
    return _pairwise_alignment("local", seq_a, seq_b, matrix_name, gap_open, gap_extend)


def _require_executable(name: str) -> str:
    exe = shutil.which(name)
    if not exe:
        raise FileNotFoundError(f"Required executable not found in PATH: {name}")
    return exe


def run_msa(records: Sequence["SeqRecord"], method: str = "muscle") -> "MultipleSeqAlignment":
    """Run MSA with an external tool and parse resulting alignment."""
    if len(records) < 2:
        raise ValueError("MSA requires at least two sequences")

    method = method.lower()
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        input_fasta = tmpdir_path / "input.fasta"
        output_fasta = tmpdir_path / "aligned.fasta"
        from Bio import SeqIO
        from Bio.Align import MultipleSeqAlignment

        SeqIO.write(records, input_fasta, "fasta")

        if method == "muscle":
            exe = _require_executable("muscle")
            cmd = [exe, "-align", str(input_fasta), "-output", str(output_fasta)]
        elif method == "mafft":
            exe = _require_executable("mafft")
            cmd = [exe, "--auto", str(input_fasta)]
        elif method == "clustalw":
            exe = _require_executable("clustalw")
            cmd = [exe, f"-INFILE={input_fasta}", f"-OUTFILE={output_fasta}", "-OUTPUT=FASTA"]
        else:
            raise ValueError("method must be one of: muscle, mafft, clustalw")

        if method == "mafft":
            with output_fasta.open("w", encoding="utf-8") as fh:
                subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.PIPE, text=True)
        else:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        aligned_records = list(SeqIO.parse(output_fasta, "fasta"))
        if not aligned_records:
            raise RuntimeError(f"MSA tool '{method}' returned no aligned sequences")
        return MultipleSeqAlignment(aligned_records)


def benchmark_msa(records: Sequence["SeqRecord"], method: str = "muscle") -> BenchmarkResult:
    tracemalloc.start()
    start = time.perf_counter()
    _ = run_msa(records, method=method)
    elapsed = time.perf_counter() - start
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return BenchmarkResult(
        method=method,
        elapsed_seconds=elapsed,
        peak_memory_kb=peak / 1024,
        n_sequences=len(records),
    )


def consensus_sequence(alignment: "MultipleSeqAlignment", threshold: float = 0.7) -> str:
    """Generate a consensus sequence using a simple per-column threshold."""
    if not alignment:
        raise ValueError("Alignment cannot be empty")

    n_rows = len(alignment)
    n_cols = alignment.get_alignment_length()
    consensus_chars: list[str] = []

    for col in range(n_cols):
        counts: dict[str, int] = {}
        for row in range(n_rows):
            residue = alignment[row, col]
            if residue == "-":
                continue
            counts[residue] = counts.get(residue, 0) + 1

        if not counts:
            consensus_chars.append("-")
            continue

        residue, count = max(counts.items(), key=lambda item: item[1])
        if (count / n_rows) >= threshold:
            consensus_chars.append(residue)
        else:
            consensus_chars.append("X")

    return "".join(consensus_chars)


def build_psiblast_command(query_fasta: str, db_name: str, iterations: int = 3) -> list[str]:
    """Return a PSI-BLAST command suitable for execution in shell."""
    if iterations < 1:
        raise ValueError("iterations must be >= 1")

    return [
        "psiblast",
        "-query",
        query_fasta,
        "-db",
        db_name,
        "-num_iterations",
        str(iterations),
        "-outfmt",
        "5",
    ]


def structural_alignment_rmsd(pdb_a: str, pdb_b: str, chain_id: str = "A") -> float:
    """Compute RMSD between alpha carbons of two protein structures."""
    from Bio.PDB import PDBParser, Superimposer

    parser = PDBParser(QUIET=True)
    struct_a = parser.get_structure("A", pdb_a)
    struct_b = parser.get_structure("B", pdb_b)

    atoms_a = [res["CA"] for res in struct_a[0][chain_id] if "CA" in res]
    atoms_b = [res["CA"] for res in struct_b[0][chain_id] if "CA" in res]

    n = min(len(atoms_a), len(atoms_b))
    if n == 0:
        raise ValueError("No alpha carbon atoms found for alignment")

    sup = Superimposer()
    sup.set_atoms(atoms_a[:n], atoms_b[:n])
    return float(sup.rms)


def to_records(ids_and_sequences: Iterable[tuple[str, str]]) -> list["SeqRecord"]:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    return [SeqRecord(Seq(preprocess_sequence(seq)), id=seq_id, description="") for seq_id, seq in ids_and_sequences]
