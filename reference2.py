# Install semanticscholar library if not already installed
#!pip install semanticscholar

# Import necessary libraries
import os
from semanticscholar import SemanticScholar
import pandas as pd
from tqdm import tqdm  # Import tqdm

# Initialize SemanticScholar client
sch = SemanticScholar()

# Function to search papers on Semantic Scholar and retrieve details
def search_and_fetch_details(query, limit=50):
    search_results = sch.search_paper(query, limit=limit)
    paper_details = []
   
    paper_count = 0
    # Using tqdm to create a progress bar
    for paper in tqdm(search_results, total=limit, desc="Fetching Papers"):
        try:
            paper_id = paper.paperId
            paper_info = sch.get_paper(paper_id)
            authors = ", ".join([author.name for author in paper_info.authors])
            citations = len(paper_info.citations) if paper_info.citations else 0

            paper_details.append({
                'title': paper_info.title,
                'abstract': paper_info.abstract,
                'authors': authors,
                'citations': citations
            })
           
        except Exception as e:
            print(f"Exception Occurred while retrieving paper with id: {paper_id}. Error: {e}")
           
        paper_count += 1
        if(paper_count >= limit):
            break

    return paper_details

# Search query
query = "Low BDDE Fillers"

# Fetch details from Semantic Scholar
papers = search_and_fetch_details(query)

# Convert the results to a DataFrame
df_papers = pd.DataFrame(papers)

# Define the path for the output CSV file on the desktop
desktop_path = os.path.join(os.path.expanduser('~'), 'Documents')
output_csv_file_path = os.path.join(desktop_path, 'semantischolartest.csv')

# Save the DataFrame to a CSV file
df_papers.to_csv(output_csv_file_path, index=False)

print(f"Data saved to {output_csv_file_path}")
