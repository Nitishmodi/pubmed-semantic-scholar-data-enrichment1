# Project Title: Automated Paper Search and Data Enrichment

## Description:
This project automates the process of searching for academic papers related to specific queries using both PubMed and Semantic Scholar APIs. It retrieves details such as title, abstract, authors, source, and URL for each paper, and saves the results into a CSV file. The project utilizes concurrent processing to improve efficiency and handles errors gracefully.

## Code Overview:
- The code combines functionalities from two references:
  - **Reference 1**: Utilizes PubMed and Semantic Scholar APIs to search for papers, fetch details, and save results into a CSV file.
  - **Reference 2**: Focuses on using the Semantic Scholar API to search for papers, fetch details, and save results into a CSV file.

## Changes Made:
- Integrated functionalities from both references to search papers from both PubMed and Semantic Scholar APIs.
- Modified the code to handle concurrent processing using ThreadPoolExecutor.
- Adjusted the code to retrieve details such as title, abstract, authors, source, and URL for each paper.
- Saved the combined results into a CSV file named 'combined_results.csv'.
- Updated the code to handle errors gracefully and provide informative error messages.
- Cleaned up the code and added comments for clarity.

## Instructions for Use:
1. Ensure Python is installed on your system.
2. Install necessary libraries using `pip install -r requirements.txt`.
3. Provide the queries in a CSV file named 'queries.csv'.
4. Run the main Python script.
5. Once completed, the combined results will be saved to a CSV file named 'combined_results.csv' in the Documents directory.

## Dependencies:
- Python 3.x
- semanticscholar
- tqdm
- pandas
- Bio (for PubMed API)

## Author:
Nitish Modi
