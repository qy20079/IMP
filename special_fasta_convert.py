from Bio import SeqIO
import random

# Function to modify the heading based on score and dataset type
def modify_heading(original_heading, score, dataset_type):
    return f"{original_heading}|{score}|{dataset_type}"

# Function to process existing FASTA file and modify headings
def modify_fasta_headings(input_file, output_file):
    with open(output_file, "w") as output_handle:
        records = list(SeqIO.parse(input_file, "fasta"))
        total_sequences = len(records)
        half_sequences = total_sequences // 2
        for sequence_number, record in enumerate(records, start=1):
            dataset_type = "training"
            score = random.randint(0, 1)  # Assigning random score of 0 or 1
            modified_heading = modify_heading(record.id, score, dataset_type)
            record.id = modified_heading
            record.description = ""
            SeqIO.write(record, output_handle, "fasta")

# Input and output file paths
input_file = "unique_seq_plus_Met.fasta"
output_file = "special_fasta_unique_seq_plus_Met.fasta"

# Modify headings and write to output file
modify_fasta_headings(input_file, output_file)