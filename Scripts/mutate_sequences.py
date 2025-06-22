from Bio import SeqIO
from Bio.Seq import Seq
import random

def mutate_sequence(seq, sub_rate=0.1, ins_rate=0.05, del_rate=0.05):
    """Introduce substitutions (10%), insertions (5%), and deletions (5%)."""
    seq = list(seq)
    new_seq = []
    i = 0
    while i < len(seq):
        # Substitutions
        if random.random() < sub_rate:
            new_seq.append(random.choice('ATCG'))
            i += 1
        # Insertions
        elif random.random() < ins_rate:
            new_seq.append(random.choice('ATCG'))  # Insert random base
        # Deletions
        elif random.random() < del_rate:
            i += 1  # Skip the current base
        else:
            new_seq.append(seq[i])
            i += 1
    return ''.join(new_seq)

# Process sequences
input_file = "last_5k_sequences.fasta"
output_file = "synthetic_non_resistant.fasta"

with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
    for record in SeqIO.parse(f_in, "fasta"):
        mutated_seq = mutate_sequence(str(record.seq))
        record.seq = Seq(mutated_seq)
        SeqIO.write(record, f_out, "fasta")
