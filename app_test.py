import streamlit as st
from streamlit_option_menu import option_menu
import xml.etree.ElementTree as ET
import xml.etree.ElementTree as ET
from src.search_engine import parse_xml
from src.search_engine import search_and_highlight
from src.search_engine import *
import re
from Bio import Entrez
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Create an option menu for the main menu in the sidebar
st.set_page_config(page_title="Search PubMed Articles", page_icon="image/logo_csie2.png")
# st.image("image/title_search.png")
st.sidebar.image("image/logo_NCKU.jpeg", use_column_width=True)
with st.sidebar:
    selected = option_menu("Main Menu", ["Search Engine", "Upload File", "Download PubMed", "Frequency Analyst"],
                           icons=['search-heart-fill', 'cloud-upload-fill', 'cloud-arrow-down-fill', 'lightbulb-fill'], menu_icon="bars", default_index=0)
# Based on the selected option, you can display different content in your web application
if selected == "Search Engine":
    st.image("image/search_title.png")
    # Sidebar

    st.sidebar.title("File Management")
    uploaded_files = st.sidebar.file_uploader("Upload PubMed XML Files", type=["xml"], accept_multiple_files=True)

    # Page

    # # Page
    # st.title("PubMed Document Retrieval")
    # st.subheader("Search PubMed Articles")

    # Initialize data list
    data = []

    # Load uploaded files
    for uploaded_file in uploaded_files:
        data += parse_xml(uploaded_file)

    # Search
    keyword = st.text_input("Enter keyword to search for:")
    case_sensitive = st.toggle("Case Sensitive Search", value=True)
    matching_articles = []

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

elif selected == "Upload File":
    st.image("image/upload_file_title.png")

elif selected == "Download PubMed":
    st.image("image/download_title.png")
    keyword_search = st.text_input("Enter keyword to download PubMed for:")

    # Add an input field for "Number of Documents" with validation
    number_of_document = st.number_input("Number of Documents", min_value=1, step=1, format="%d")

    if st.button("Download"):
        folder_path = f"dataset/{keyword_search}"
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        else:
            st.warning(f"Folder '{folder_path}' already exists.", icon="⚠️")

        # Search PubMed database
        handle = Entrez.esearch(db="pubmed", term=keyword_search, retmax=number_of_document)
        record = Entrez.read(handle)
        handle.close()

        # Download and save XML files
        progress_text = "Please wait! Downloading ..."
        my_bar = st.progress(0, text=progress_text)
        number_of_document_real = len(record["IdList"])
        process = 0
        for i, pubmed_id in enumerate(tqdm(record["IdList"], desc="Processing documents")):
            fetch_handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="xml", retmode="xml")
            xml_data = fetch_handle.read()
            fetch_handle.close()

            # Construct XML file's name, typically using PubMed ID
            xml_file_name = os.path.join(folder_path, f"{pubmed_id}.xml")

            # Save XML data to a file, using binary mode 'wb' to write byte data
            with open(xml_file_name, "wb") as xml_file:
                xml_file.write(xml_data)
            process += int(100/number_of_document_real)
            my_bar.progress( process  , text=progress_text)
        my_bar.empty()
        st.success(' Download success!', icon="✅")
        st.balloons()

    path_keywords = os.listdir("dataset")
    if len(path_keywords) > 0: 
        number_of_files = []
        for path_keyword in path_keywords:
            folder_path = f"dataset/{path_keyword}" 
            file_count = len(os.listdir(folder_path))
            #print(folder_path)
            print(file_count)
            number_of_files.append(file_count) 


        # Create a list of dictionaries to store the data
        st.markdown(f'<p style="text-align:center; color:red;">The number of articles containing the keyword in the dataset</p>', unsafe_allow_html=True)
        table_data = [{"Keyword": path_keywords[i], 
                    "Number of articles": number_of_files[i]} for i in range(len(path_keywords))]
        st.table(table_data)

elif selected == "Frequency Analyst":
    st.image("image/analyst_title.png")
    keyword_search = ''
    keyword_search = st.text_input("Enter keyword :")
    edit_distance = st.toggle("Edit distance", value=False)
    path_keywords = os.listdir("dataset")
    if keyword_search in path_keywords: 
        st.info(f"Your keyword is {keyword_search}", icon="ℹ️")
    elif edit_distance and len(keyword_search) > 0: 
        suggestions = find_closest_keywords(keyword_search, path_keywords, num_suggestions = 10  )

        # Create a bar chart using Seaborn
        keywords, probabilities = zip(*suggestions)
        data = {"Keywords": keywords, "Probability": probabilities}
        df = pd.DataFrame(data)

        plt.figure(figsize=(8, 4))
        sns.set(style="whitegrid")  # Set the style to have a white grid
        ax = sns.barplot(x="Keywords", y="Probability", data=df)
        plt.xticks(rotation=45)
        plt.title("Keyword Probabilities with Edit Distance")

        # Add percentages on each bar
        for p in ax.patches:
            ax.annotate(f'{p.get_height()*100:.2f}%', (p.get_x() + p.get_width() / 2., p.get_height()),
                        ha='center', va='center', fontsize=10, color='black', xytext=(0, 5),
                        textcoords='offset points')

        st.pyplot(plt)
        # st.write(f"Closest keywords to '{keyword_search}' is {keywords[0]}") 
        if probabilities[0] > 0.7: 
            st.info(f"Closest keywords to '{keyword_search}': {keywords[0]}",icon="ℹ️")
        else:
            st.warning("No found the keyword", icon="⚠️")
            keyword_search = '' 
    elif len(keyword_search) > 0: 
        st.warning("No found the keyword", icon="⚠️")
        keyword_search = ''
    if len(keyword_search) > 0: 
        st.sidebar.title("Setting")
        top_of_word = st.sidebar.number_input("Top of words", min_value=1, step=1, format="%d")

        if st.sidebar.button('Start analyzing ...'):
            st.write("Start analyzing ...")








