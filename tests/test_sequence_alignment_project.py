import pytest

from src.sequence_alignment_project import (
    build_psiblast_command,
    consensus_sequence,
    pairwise_global_alignment,
    pairwise_local_alignment,
    preprocess_sequence,
    to_records,
)


def test_preprocess_sequence_removes_invalid_characters() -> None:
    assert preprocess_sequence("mk-ta*yi??") == "MKTAYI"


def test_build_psiblast_command() -> None:
    cmd = build_psiblast_command("query.fasta", "nr", iterations=5)
    assert cmd[0] == "psiblast"
    assert "-num_iterations" in cmd
    assert "5" in cmd


def test_build_psiblast_command_requires_positive_iterations() -> None:
    with pytest.raises(ValueError):
        build_psiblast_command("query.fasta", "nr", iterations=0)


def test_pairwise_global_alignment_metrics() -> None:
    pytest.importorskip("Bio")

    result = pairwise_global_alignment("MKTAYI", "MKTQYI")
    assert result.score > 0
    assert len(result.aligned_seq_a) == len(result.aligned_seq_b)
    assert 0 <= result.identity <= 1
    assert 0 <= result.similarity <= 1


def test_pairwise_local_alignment_metrics() -> None:
    pytest.importorskip("Bio")

    result = pairwise_local_alignment("AAAAAMKTAYI", "MKTQYI")
    assert result.score > 0
    assert len(result.aligned_seq_a) == len(result.aligned_seq_b)
    assert result.identity >= 0


def test_pairwise_alignment_rejects_invalid_sequences() -> None:
    with pytest.raises(ValueError):
        pairwise_global_alignment("---***", "???")


def test_consensus_sequence() -> None:
    pytest.importorskip("Bio")
    from Bio.Align import MultipleSeqAlignment

    records = to_records([
        ("s1", "MKTAYI"),
        ("s2", "MKTQYI"),
        ("s3", "MKTAYI"),
    ])
    alignment = MultipleSeqAlignment(records)
    consensus = consensus_sequence(alignment, threshold=0.6)
    assert consensus.startswith("MKT")
