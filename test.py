from Bio import Entrez
from Bio import Medline
import os

# Set your email address for NCBI
Entrez.email = "ben@iir.csie.ncku.edu.tw"

# Your search term
search_term = "lymphoma"

# Number of articles to retrieve
num_articles = 5

# Perform the search
handle = Entrez.esearch(db="pubmed", term=search_term, retmax=num_articles)
record = Entrez.read(handle)
print(record)
handle.close()

# Fetch and save the articles as XML files
for article_id in record["IdList"]:
    handle = Entrez.efetch(db="pubmed", id=article_id, rettype="medline", retmode="xml")
    records = Medline.parse(handle)
    for record in records:
        if 'PMID' in record:
            filename = f"{record['PMID']}.xml"
            with open(filename, "w") as file:
                file.write(handle.read())
                print(handle.read())
    handle.close()

