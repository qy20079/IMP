import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Function to generate random protein sequence of a given length with a limited set of amino acids
def generate_random_sequence(length):
    bases = "ACD"  # Limited set of amino acids
    return ''.join(random.choice(bases) for _ in range(length))

# Function to determine dataset type (training or testing) based on sequence number
def get_dataset_type(sequence_number, total_sequences):
    return "training" if sequence_number <= total_sequences / 2 else "testing"

# Function to determine score (0 or 1) based on sequence number
def get_score(sequence_number, total_sequences):
    return 0 if sequence_number <= total_sequences / 4 or (sequence_number > total_sequences / 2 and sequence_number <= total_sequences * 3 / 4) else 1

# Load existing sequences
existing_sequences = list(SeqIO.parse("unique_seq_plus_Met.fasta", "fasta"))

# Calculate total number of existing sequences
total_existing_sequences = len(existing_sequences)

# total number of new sequences to add
total_new_sequences = 500

# Total number of sequences after adding new ones
total_sequences = total_existing_sequences + total_new_sequences

# Shuffle existing sequences to assign random scores and dataset types
random.shuffle(existing_sequences)

# Determine dataset type and score for existing sequences
for i, seq_record in enumerate(existing_sequences, start=1):
    dataset_type = get_dataset_type(i, total_sequences)
    score = get_score(i, total_sequences)
    heading = f"{seq_record.id}|{score}|{dataset_type}"
    seq_record.id = heading

# Generate and add new sequences
for i in range(total_existing_sequences + 1, total_sequences + 1):
    dataset_type = "training" if i <= total_sequences / 2 else "testing"
    score = 1  # All new sequences are given a score of 1
    seq = generate_random_sequence(50)  # Generating a random sequence of length 50
    heading = f"New_Sequence_{i}|{score}|{dataset_type}"
    seq_record = SeqRecord(Seq(seq), id=heading, description="")
    existing_sequences.append(seq_record)

# Write all sequences (existing + new) to a new FASTA file
with open("combined_random_existing_sequences.fasta", "w") as output_handle:
    SeqIO.write(existing_sequences, output_handle, "fasta")
