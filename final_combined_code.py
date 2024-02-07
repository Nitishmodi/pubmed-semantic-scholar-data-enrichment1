import csv
import pandas as pd
import random
import time
import os
from Bio import Entrez
from semanticscholar import SemanticScholar
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# Set the email here to inform NCBI who is making the request
Entrez.email = "vihangsharma@gmail.com"

# Initialize SemanticScholar client
sch = SemanticScholar()

def search_pubmed(query):
    """Search PubMed for the given query and return the list of articles."""
    try:
        handle = Entrez.esearch(db='pubmed', sort='relevance', retmax=20, retmode='xml', term=query)
        results = Entrez.read(handle)
        handle.close()
        id_list = results['IdList']
        if not id_list:
            return []
        handle = Entrez.efetch(db='pubmed', retmode='xml', id=','.join(id_list))
        details = Entrez.read(handle)
        handle.close()
        pubmed_results = []
        for paper in details['PubmedArticle']:
            article = paper['MedlineCitation']['Article']
            title = article.get('ArticleTitle', '')
            abstract = article.get('Abstract', {}).get('AbstractText', '')
            authors_list = article.get('AuthorList', [])
            authors = ', '.join([author.get('LastName', '') + ' ' + author.get('ForeName', '') for author in authors_list if 'LastName' in author and 'ForeName' in author])
            pmid = paper['MedlineCitation']['PMID']
            url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            pubmed_results.append({'title': title, 'abstract': abstract, 'authors': authors, 'source': 'PubMed', 'url': url})
        return pubmed_results
    except Exception as e:
        print(f"Error in PubMed search: {e}")
        return []

def search_semantic_scholar(query, limit=5):
    """Search Semantic Scholar for the given query and fetch details of articles."""
    try:
        search_results = sch.search_paper(query, limit=limit)
        semantic_results = []
        for paper in tqdm(search_results, total=limit, desc="Fetching Papers"):
            try:
                paper_id = paper.paperId
                paper_info = sch.get_paper(paper_id)
                authors = ", ".join([author.name for author in paper_info.authors])
                semantic_results.append({'title': paper_info.title, 'abstract': paper_info.abstract, 'authors': authors, 'source': 'Semantic Scholar', 'url': paper_info.url})
            except Exception as e:
                print(f"Exception Occurred while retrieving paper with id: {paper_id}. Error: {e}")
        return semantic_results
    except Exception as e:
        print(f"Error in Semantic Scholar search: {e}")
        return []

def process_query(query):
    """Process a single query through both PubMed and Semantic Scholar."""
    combined_results = search_pubmed(query)
    combined_results.extend(search_semantic_scholar(query, limit=1))
    # Randomized sleep between 0.3 and 0.7 seconds to respect rate limits
    time.sleep(random.uniform(0.3, 0.7))
    return combined_results

def main():
    # Load queries from 'queries.csv'
    try:
        queries = []
        with open('queries.csv', 'r', newline='') as infile:
            reader = csv.reader(infile)
            queries = [row[0] for row in reader]
    except Exception as e:
        print(f"Error reading queries: {e}")
        return

    # Process queries concurrently, up to 4 at a time
    combined_results = []
    with ThreadPoolExecutor(max_workers=4) as executor:
        future_to_query = {executor.submit(process_query, query): query for query in queries}
        for future in as_completed(future_to_query):
            combined_results.extend(future.result())

    # Convert combined results to a DataFrame and save to CSV
    df = pd.DataFrame(combined_results)
    
    # Define the path for the output CSV file
    output_csv_file_path = os.path.join(os.path.expanduser('~'), 'Documents', 'combined_results.csv')
    
    try:
        df.to_csv(output_csv_file_path, index=False)
        print(f"CSV file generated successfully at: {output_csv_file_path}")
    except Exception as e:
        print(f"Error saving to CSV: {e}")

if __name__ == "__main__":
    main()
