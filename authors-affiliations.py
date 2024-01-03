
import csv
from Bio import Entrez

def get_paper_info(pubmed_id):
    Entrez.email = "your_email@example.com"  # Set your email address

    # Fetch the paper information from PubMed
    handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="medline", retmode="text")
    record = Entrez.read(handle)

    # Extract authors and affiliations
    authors = record[0]['AU']
    affiliations = [author_info.get('AD', '') for author_info in record[0].get('FAU', [])]

    return authors, affiliations

def write_to_csv(pubmed_id, authors, affiliations):
    filename = f"pubmed_{pubmed_id}_authors_affiliations.csv"

    with open(filename, mode='w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write header
        writer.writerow(['Author', 'Affiliation'])

        # Write data
        for author, affiliation in zip(authors, affiliations):
            writer.writerow([author, affiliation])

    print(f"Data written to {filename}")

def main():
    # Get PubMed ID from user input
    pubmed_id = input("Enter PubMed ID: ")

    # Get and print paper information
    authors, affiliations = get_paper_info(pubmed_id)

    # Write to CSV
    write_to_csv(pubmed_id, authors, affiliations)

if __name__ == "__main__":
    main()

