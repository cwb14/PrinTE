import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Function to insert the tag into sequences and update headers
def insert_tag_and_update_header(input_fasta, output_fasta, tag_sequence='CCCCCCCCCCCCCCCCCC'):
    unassigned_count = 0  # Counter for unassigned sequences
    with open(output_fasta, 'w') as outfile:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            seq = str(record.seq)
            length = len(seq)
            midpoint = length // 2  # Integer division for midpoint

            # Insert the tag at the midpoint
            new_seq = seq[:midpoint] + tag_sequence + seq[midpoint:]

            # Update the header based on suffix
            if '_intact' in record.id or '_solo' in record.id:
                new_id = record.id  # No changes needed
            else:
                new_id = f"{record.id}_intact"
                unassigned_count += 1

            # Create a new SeqRecord with updated sequence and ID
            updated_record = SeqRecord(
                id=new_id,
                description="",  # Clear the description
                seq=new_seq
            )

            # Write the updated record to the output file
            SeqIO.write(updated_record, outfile, 'fasta')

    print(f"Tag inserted and headers updated. Output written to '{output_fasta}'.")
    print(f"{unassigned_count} unassigned sequences were assumed to be intact and tagged with '_intact'.")

# Main function for argument parsing
def main():
    parser = argparse.ArgumentParser(description="Insert a tag into the midpoint of each sequence in a FASTA file and update headers.")
    parser.add_argument('-i', '--input', required=True, help="Input FASTA file")
    parser.add_argument('-o', '--output', required=True, help="Output FASTA file")
    args = parser.parse_args()

    # Call the tag insertion and header update function
    insert_tag_and_update_header(args.input, args.output)

if __name__ == "__main__":
    main()
