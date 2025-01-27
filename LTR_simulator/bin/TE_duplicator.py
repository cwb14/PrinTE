from Bio import SeqIO

# Input and output file names
input_fasta = "curated_LTRs_young_proof.fa"
output_fasta = "curated_LTRs_young_proof_dup.fa"

# Number of duplicates for each sequence
num_duplicates = 50

# Read the input FASTA file and duplicate sequences
with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        for i in range(1, num_duplicates + 1):
            # Modify the header to include the duplicate index
            new_id = f"{record.id}_{i}"
            new_description = f"{record.description}_{i}"
            new_record = record[:]
            new_record.id = new_id
            new_record.description = new_description
            # Write the modified sequence to the output file
            SeqIO.write(new_record, outfile, "fasta")

print(f"Duplicated sequences written to {output_fasta}")
