from Bio import Entrez

Entrez.email = "prasanth.thuthika00@gmail.com.com"  

# Search for E. coli CDS sequences (RefSeq)
search_term = "Escherichia coli[orgn] AND CDS[feature] AND refseq[filter]"
search_handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=10000)
search_results = Entrez.read(search_handle)
search_handle.close()

# Fetch sequences
id_list = search_results["IdList"]
fetch_handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
fasta_data = fetch_handle.read()
fetch_handle.close()

# Save to file
with open("ecoli_cds.fasta", "w") as f:
    f.write(fasta_data)
