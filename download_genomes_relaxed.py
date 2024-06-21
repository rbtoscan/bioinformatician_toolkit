
"""
Relaxed Genome Downloader

This script reads a list of pathogen names from a specified file and attempts to download their genomes
from the NCBI database. It first searches for complete genomes. If no complete genome is found, it will
download the best available genome assembly (complete or not). The script prints messages indicating 
whether a complete genome or the best available genome was downloaded, along with the Genome ID.

Usage:
    python download_genomes.py <organism_list_file>

Example:
     python download_genomes.py list_genomes.txt

"""


from Bio import Entrez, SeqIO
import os
import sys

# Email to identify yourself to NCBI
Entrez.email = "your_email@example.com"

# Function to read organisms from a file
def read_organism_list(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file]

# Output directory
output_dir = "downloaded_genomes"
os.makedirs(output_dir, exist_ok=True)

def download_genome(organism):
    # Try to find complete genome
    search_handle = Entrez.esearch(db="nucleotide", term=f"{organism}[ORGN] AND complete genome[TI]", retmax=1)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    genome_type = "complete genome"
    if not search_results["IdList"]:
        # If no complete genome found, search for the best available
        search_handle = Entrez.esearch(db="nucleotide", term=f"{organism}[ORGN]", retmax=1)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        genome_type = "best available genome assembly (complete or not)"
    
    if search_results["IdList"]:
        genome_id = search_results["IdList"][0]
        fetch_handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")
        genome_data = fetch_handle.read()
        fetch_handle.close()
        
        # Save the genome data to a file
        file_path = os.path.join(output_dir, f"{organism.replace(' ', '_')}.fasta")
        with open(file_path, "w") as genome_file:
            genome_file.write(genome_data)
        print(f"Downloaded {genome_type} for {organism} (Genome ID: {genome_id})")
    else:
        print(f"No genome found for {organism}")

# Read organisms from a file
organism_list_file = sys.argv[1]
organisms = read_organism_list(organism_list_file)

# Download genomes for all organisms in the list
for organism in organisms:
    download_genome(organism)
