from Bio import SeqIO

def extract_fasta_sequence(input_fasta, start, end, output_fasta):
    """
    Extracts a subsequence from a FASTA file based on the start and end positions.

    Parameters:
    input_fasta (str): Path to the input FASTA file.
    start (int): Start position (1-based indexing).
    end (int): End position (1-based indexing).
    output_fasta (str): Path to the output FASTA file.
    """
    # Read the FASTA file
    with open(input_fasta, "r") as input_handle:
        # Parse the FASTA file
        for record in SeqIO.parse(input_handle, "fasta"):
            # Extract the subsequence
            subsequence = record.seq[start-1:end]  # Python indexing starts at 0
            # Create a newrecord with the extracted subsequence
            new_record = record[start-1:end]
            new_record.id = f"{record.id}_extract_{start}_{end}"
            new_record.description = f"Extracted nucleotides from {start} to {end}"

            # Write the subsequence to the output file
            with open(output_fasta, "w") as output_handle:
                SeqIO.write(new_record, output_handle, "fasta")

    print(f"Nucleotides from {start} to {end} have been extracted and saved to {output_fasta}.")

# Example usage:
# extract_fasta_sequence("MhaplaVW9.mitogenome_linear.fa", 11535, 11580, "trnR_hap.fasta")
# extract_fasta_sequence("MhaplaVW9.mitogenome_linear.fa", 11582, 11633, "trnA_hap.fasta")
# extract_fasta_sequence("MhaplaVW9.mitogenome_linear.fa", 11634, 11687, "trnQ_hap.fasta")


extract_fasta_sequence("MhaplaVW9.mitogenome_linear.fa",4214 ,4415, "betweenrrnaL.fasta")