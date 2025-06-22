from Bio import Entrez

Entrez.email = "prasanth.thuthika00@gmail.com" 

search_term = "Bacillus subtilis[orgn] AND CDS[feature] AND refseq[filter]"
search_handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=10000)
search_results = Entrez.read(search_handle)
search_handle.close()

id_list = search_results["IdList"]
fetch_handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
fasta_data = fetch_handle.read()
fetch_handle.close()


with open("bacillus_cds.fasta", "w") as f:
    f.write(fasta_data)
