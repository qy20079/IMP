from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random

# Function to generate a random nucleotide sequence of a given length
def generate_random_sequence(length):
    nucleotides = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(nucleotides) for _ in range(length))

# Function to generate multiple random sequences and append them to the original sequence
def combine_sequences(input_fasta, output_fasta, num_random_sequences, sequence_length):
    # Read the existing sequence from the input FASTA file
    records = list(SeqIO.parse(input_fasta, "fasta"))
    existing_sequence = str(records[0].seq)
    
    # Generate the specified number of random nucleotide sequences and combine them
    random_sequences = ''.join(generate_random_sequence(sequence_length) for _ in range(num_random_sequences))
    
    # Combine the existing sequence with the generated random sequences
    combined_sequence = existing_sequence + random_sequences
    
    # Create a new SeqRecord object for the combined sequence
    combined_record = SeqRecord(Seq(combined_sequence), id="combined_sequence", description="Original + 300 Random Sequences")
    
    # Write the combined sequence to the output FASTA file
    SeqIO.write(combined_record, output_fasta, "fasta")
    print(f"Combined sequence saved to {output_fasta}")

# Specify input and output file paths
input_fasta = "existing_sequence.fasta"  # Path to the existing FASTA file
output_fasta = "biased_nt.fasta"  # Output file for the combined sequence

# Number of random sequences and their length
num_random_sequences = 300  # Number of random sequences
sequence_length = 100       # Length of each random sequence

# Run the sequence combination
combine_sequences(input_fasta, output_fasta, num_random_sequences, sequence_length)
ombine_sequences(input_fasta, output_fasta)
