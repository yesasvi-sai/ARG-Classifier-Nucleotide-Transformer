import csv
from Bio import SeqIO

def process_fasta_to_csv(fasta_file, label):
    csv_file = fasta_file.replace('.fasta', '.csv')
    
    # Extract field names from the header (using the first sequence as reference)
    first_record = next(SeqIO.parse(fasta_file, "fasta"))
    header_fields = [field.split(':')[0] for field in first_record.description.split('|')]
    
    # Define CSV columns (header fields + sequence + label)
    csv_columns = header_fields + ["sequence", "label"]
    
    # Write to CSV
    with open(csv_file, 'w', newline='') as csv_out:
        writer = csv.DictWriter(csv_out, fieldnames=csv_columns)
        writer.writeheader()
        
        for record in SeqIO.parse(fasta_file, "fasta"):
            row_data = {}
            # Parse header fields
            for field in record.description.split('|'):
                key, value = field.split(':', 1)  # Split into key-value pair
                row_data[key] = value
            # Add sequence and label
            row_data["sequence"] = str(record.seq)
            row_data["label"] = label
            writer.writerow(row_data)
    
    print(f"Created {csv_file} with {len(list(SeqIO.parse(fasta_file, 'fasta')))} sequences")

# Process resistant sequences (label=1)
process_fasta_to_csv("resistant_sequences.fasta", label=1)

# Process synthetic non-resistant sequences (label=0)
process_fasta_to_csv("synthetic_non_resistant.fasta", label=0)
