from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random

# Generate random nucleotide sequence of a given length
def generate_random_sequence(length):
    bases = "ATCG"  # Nucleotides
    return ''.join(random.choice(bases) for _ in range(length))

# Generate random set of nucleotide sequences with special fasta headings
def generate_sequences_with_headings(num_sequences, sequence_length):
    sequences = []
    for i in range(1, num_sequences + 1):
        # Generate random nucleotide sequence
        seq = generate_random_sequence(sequence_length)
        # Determine score (0 or 1)
        score = i % 2
        # Determine dataset type (training or testing)
        dataset_type = "training" if i <= num_sequences / 2 else "testing"
        # Create heading
        heading = f"|{score}|{dataset_type}"
        # Create SeqRecord
        seq_record = SeqRecord(Seq(seq), id=f"Random_Sequence_{i}", description=heading)
        sequences.append(seq_record)
    return sequences

# Parameters
num_sequences = 1000
sequence_length = 100

# Generate sequences with headings
sequences = generate_sequences_with_headings(num_sequences, sequence_length)

# Save sequences to a FASTA file
SeqIO.write(sequences, "random_nucleotide_sequences.fasta", "fasta")

# Print sequences
for seq in sequences:
    print(seq.format("fasta"))
