import pandas as pd
import io
import json
from Bio import Entrez
import streamlit as st
from dotenv import load_dotenv
import os
import re

# Load environment variables from the .env file
load_dotenv()

# Set email and API key
Entrez.email = 'saikatsinha21@gmail.com'
Entrez.api_key = os.getenv('ENTREZ_API_KEY')

# Streamlit Web App
st.title("PubMed Article Search")

# Input fields
authors_input = st.text_input("Enter authors (comma-separated)", "")
topics_input = st.text_input("Enter topics (comma-separated)", "")
start_date = st.date_input("Start date", value=pd.to_datetime("2012-03-01"))
end_date = st.date_input("End date", value=pd.to_datetime("2022-12-31"))

# Process user inputs
authors = [author.strip() for author in authors_input.split(",") if author.strip()]
topics = [topic.strip() for topic in topics_input.split(",") if topic.strip()]
date_range = f'("{start_date.strftime("%Y/%m/%d")}"[Date - Create] : "{end_date.strftime("%Y/%m/%d")}"[Date - Create])'

# Build the query dynamically
queries = []
if authors:
    author_queries = ['{}[Author]'.format(author) for author in authors]
    queries.append('(' + ' OR '.join(author_queries) + ')')

if topics:
    topic_queries = ['{}[Title/Abstract]'.format(topic) for topic in topics]
    queries.append('(' + ' OR '.join(topic_queries) + ')')

full_query = ' AND '.join(queries) + ' AND ' + date_range

# Search and Fetch Articles
if st.button("Search PubMed"):
    handle = Entrez.esearch(db='pubmed', retmax=50, term=full_query)
    record = Entrez.read(handle)
    id_list = record['IdList']

    # DataFrame to store results
    df = pd.DataFrame(columns=['PubmedID', 'Title', 'Publication Date', 'Non-academic Author(s)', 
                               'Company Affiliation(s)', 'Corresponding Author Email'])

    for pmid in id_list:
        handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
        records = Entrez.read(handle)

        for record in records['PubmedArticle']:
            title = record['MedlineCitation']['Article']['ArticleTitle']

            # Extract publication date
            pub_date = record['MedlineCitation']['Article']['Journal']['JournalIssue'].get('PubDate', {})
            year = pub_date.get('Year', '')
            month = pub_date.get('Month', '')
            day = pub_date.get('Day', '')
            month_mapping = {
                "Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04", "May": "05", "Jun": "06",
                "Jul": "07", "Aug": "08", "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12"
            }
            month = month_mapping.get(month, month)
            if year:
                if month and day:
                    publication_date = f"{year}-{month}-{day}"
                elif month:
                    publication_date = f"{year}-{month}-01"
                else:
                    publication_date = f"{year}-01-01"
            else:
                publication_date = ''

            # Extract authors and affiliations
            authors = []
            non_academic_authors = []
            company_affiliations = []
            corresponding_author_email = ''
            
            if 'AuthorList' in record['MedlineCitation']['Article']:
                for author in record['MedlineCitation']['Article']['AuthorList']:
                    if 'LastName' in author and 'ForeName' in author:
                        authors.append(f"{author['LastName']} {author['ForeName']}")
                    
                    # Check for affiliations
                    if 'AffiliationInfo' in author and author['AffiliationInfo']:
                        affiliation = author['AffiliationInfo'][0]['Affiliation']
                        if any(company in affiliation.lower() for company in ['pharma', 'biotech', 'inc', 'ltd', 'corp']):
                            company_affiliations.append(affiliation)
                        else:
                            non_academic_authors.append(f"{author['LastName']} {author['ForeName']}")

                    # Check for email in the affiliation
                    if 'AffiliationInfo' in author and author['AffiliationInfo']:
                        if '@' in author['AffiliationInfo'][0]['Affiliation']:
                            corresponding_author_email = author['AffiliationInfo'][0]['Affiliation'].split()[-1]

            # Join lists to create strings
            non_academic_authors = ', '.join(set(non_academic_authors))
            company_affiliations = ', '.join(set(company_affiliations))

            new_row = pd.DataFrame({
                'PubmedID': [pmid],
                'Title': [title],
                'Publication Date': [publication_date],
                'Non-academic Author(s)': [non_academic_authors],
                'Company Affiliation(s)': [company_affiliations],
                'Corresponding Author Email': [corresponding_author_email]
            })

            df = pd.concat([df, new_row], ignore_index=True)

    # Display DataFrame
    st.write("Search Results:", df)

    # Convert DataFrame to CSV and download
    csv = df.to_csv(index=False)
    st.download_button(
        label="Download results as CSV",
        data=csv,
        file_name='PubMed_results.csv',
        mime='text/csv'
    )