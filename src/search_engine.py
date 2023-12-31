import streamlit as st
from streamlit_option_menu import option_menu
import xml.etree.ElementTree as ET
import xml.etree.ElementTree as ET
from src.utils import parse_xml
from src.utils import search_and_highlight
from src.utils import *
import re
from Bio import Entrez
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import nltk
from nltk.corpus import stopwords
from collections import Counter
from nltk.stem import PorterStemmer
from fuzzywuzzy import fuzz
from fuzzywuzzy import process

def search_engine():
    st.image("image/search_title.png")
    # Sidebar

    # st.sidebar.title("File Management")
    # uploaded_files = st.sidebar.file_uploader("Upload PubMed XML Files", type=["xml"], accept_multiple_files=True)
    # list_file = os.listdir("dataset/Covid19") 
    # print(uploaded_files)

    # Search
    keyword = st.text_input("Enter keyword to search for:")
    case_sensitive = True
    # case_sensitive = st.toggle("Case Sensitive Search", value=True)
    matching_articles = []
    # edit_distance = st.toggle("Edit distance", value=False)
    edit_distance = True

    if edit_distance and len(keyword) > 0: 
        data = []
        documents = []
        # Load uploaded files
        filtered_tokens = [] 
        list_file = os.listdir("dataset/Covid19")
        for file in list_file:
            try:
              data += parse_xml(f"dataset/Covid19/{file}")
              documents.append(parse_xml_to_string(f"dataset/Covid19/{file}"))
            except: 
              continue 

        for doc in documents:
          tokens = clean_and_tokenize(doc)
          filtered_tokens.extend([word for word in tokens if word not in stopwords.words('english')])

        len(filtered_tokens)
        unique_list = list(set(map(str.lower, filtered_tokens)))

        filtered_list = [word for word in unique_list if not word.isnumeric() and word.isalpha() and len(word) <= 10]

        # Calculate similarity scores for all keywords
        similarity_scores = [(word, fuzz.ratio(keyword, word)) for word in filtered_list]
        # Sort the list of keywords by similarity in descending order
        sorted_keywords = sorted(similarity_scores, key=lambda x: x[1], reverse=True)
        # Get the top 10 keywords
        top_10_keywords = [word for word, score in sorted_keywords[:9]]
        top_10_keywords.append(keyword)

        option = st.selectbox('Recommend the nearest word...',top_10_keywords)

        st.info(f"Your keyword: {option}", icon="ℹ️")
        keyword = option        

        # print("Top 10 keywords:")
        # for keyword in top_10_keywords:
        #     print(keyword)
          

    if st.button("Search"):
        for article in data:
            highlighted_fields = search_and_highlight(article, keyword, case_sensitive)
            if any(isinstance(value, str) and '<span style="background-color: yellow">' in value for value in highlighted_fields.values()):
                matching_articles.append(highlighted_fields)

        if not matching_articles:
            st.write("No matching articles found.")
        else:

            for idx, article in enumerate(matching_articles, start=1):
                st.markdown(f'<p style="text-align:center; color:red;">Matching Article {idx}:</p>', unsafe_allow_html=True)
                
                # Calculate and display line count for abstract
                abstract_text = article['Abstract']
                num_lines = len(abstract_text.split('\n')) if abstract_text else 0
                # st.markdown(f"**Number of Lines in Abstract**: {num_lines}", unsafe_allow_html=True)

                # Calculate and display document statistics
                try:
                    num_characters = len(article['Abstract'])
                except TypeError:
                    num_characters = 0

                try:
                    abstract_text = article['Abstract']
                    num_words = len(abstract_text.split())
                except (TypeError, AttributeError):
                    num_words = 0

                num_sentences = len(re.split(r'[.!?]', abstract_text)) if abstract_text else 0

                # Create a table for document statistics
                statistics_table = {
                    "Statistic": ["Number of Characters", "Number of Words", "Number of Sentences (EOS)"],
                    "Value": [num_characters, num_words, num_sentences]
                }
                st.table(statistics_table)

                # Display other article information
                for key, value in article.items():
                    if key in ['PMID', 'Title', 'Journal Title', 'ISSN', 'Publication Date', 'Authors', 'Keywords']:
                        # Format these fields as bold and italic
                        st.markdown(f"**_{key}_**: {value}", unsafe_allow_html=True)
                    else:
                        st.markdown(f"**{key}**: {value}", unsafe_allow_html=True)
                st.write("---")

        # Display the total number of matching articles
        st.write(f"Total Number of Matching Articles: {len(matching_articles)}")