from pathlib import Path

from src.sequence_alignment_project import (
    consensus_sequence,
    load_fasta_sequences,
    pairwise_global_alignment,
    pairwise_local_alignment,
)


def main() -> None:
    records = load_fasta_sequences(Path("data/sample_sequences.fasta"))
    seq_a = str(records[0].seq)
    seq_b = str(records[1].seq)

    global_result = pairwise_global_alignment(seq_a, seq_b)
    local_result = pairwise_local_alignment(seq_a, seq_b)

    print("Global score:", global_result.score)
    print("Global identity:", f"{global_result.identity:.2%}")
    print("Local score:", local_result.score)
    print("Local similarity:", f"{local_result.similarity:.2%}")

    print("\nGlobal alignment block:")
    print("seq1", global_result.aligned_seq_a)
    print("seq2", global_result.aligned_seq_b)

    from Bio.Align import MultipleSeqAlignment

    alignment = MultipleSeqAlignment(records)
    print("Consensus:", consensus_sequence(alignment))


if __name__ == "__main__":
    main()
