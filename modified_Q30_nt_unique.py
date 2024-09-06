from Bio import SeqIO

# Function to convert to special fasta heading 
def modify_heading(original_heading, score, dataset_type):
    return f"{original_heading}|{score}|{dataset_type}"

# Function to process existing FASTA files and modify headings
def modify_fasta_headings(input_files, output_file, scores):
    with open(output_file, "w") as output_handle:
        for input_file, score in zip(input_files, scores):
            records = list(SeqIO.parse(input_file, "fasta"))
            for record in records:
                dataset_type = "training"
                modified_heading = modify_heading(record.id, score, dataset_type)
                record.id = modified_heading
                record.description = ""
                SeqIO.write(record, output_handle, "fasta")

# Input and output file paths
input_files = ["green_Q30_nt_unique.fasta", "August_dim_Q30_nt_unique.fasta"]
scores = [1, 0]  # Scores for the input files
output_file = "dim_green_combined_Q30_nt_unique.fasta"

# Modify headings and write to output file
modify_fasta_headings(input_files, output_file, scores)
