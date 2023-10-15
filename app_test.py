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
import nltk
from nltk.corpus import stopwords
from collections import Counter
from nltk.stem import PorterStemmer

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
            keyword_search = keywords[0]
        else:
            st.warning("No found the keyword", icon="⚠️")
            keyword_search = '' 
    elif len(keyword_search) > 0: 
        st.warning("No found the keyword", icon="⚠️")
        keyword_search = ''
    if len(keyword_search) > 0: 
        st.sidebar.title("Setting")
        top_of_word = st.sidebar.number_input("Top of words", min_value=5, step=1, format="%d")
        zip_distribution  = st.sidebar.toggle("Zipf Distribution", value=False)
        if zip_distribution:
            st.write("Start analyzing ...")

            list_file = os.listdir(f"dataset/{keyword_search}") 

            documents = []
            for file in list_file: 
                documents.append(parse_xml_to_string(f"dataset/ben/{file}")) 

            filtered_tokens = []
            for doc in documents:
                tokens = clean_and_tokenize(doc)
                filtered_tokens.extend(tokens)

            # Calculate word frequencies
            word_freq = Counter(filtered_tokens)

            # Sort the words by frequency in descending order
            sorted_word_freq = dict(sorted(word_freq.items(), key=lambda item: item[1], reverse=True))

            # Display only the top 20 frequencies
            number_of_words = top_of_word
            top_words = dict(list(sorted_word_freq.items())[:number_of_words])

            # Set Seaborn style
            sns.set(style="whitegrid")

            # Plot the Zipf distribution for the top 20 words using Seaborn
            plt.figure(figsize=(10, 6))
            ax = sns.barplot(x=list(top_words.values()), y=list(top_words.keys()), palette="viridis")
            ax.set(xlabel='Frequency', ylabel='Words')
            plt.title(f'Top {number_of_words} Words | Zipf Distribution of Terms (with stopwords)')
            plt.tight_layout()

            # Display the plot
            st.pyplot(plt)

            # Remove stopwords
            remove_stopwords  = st.sidebar.toggle("Remove Stopwords", value=False)
            if remove_stopwords:
                list_file = os.listdir(f"dataset/{keyword_search}") 

                documents = []
                for file in list_file: 
                    documents.append(parse_xml_to_string(f"dataset/ben/{file}")) 

                # Remove stopwords and tokenize the text
                filtered_tokens = []
                for doc in documents:
                    tokens = clean_and_tokenize(doc)
                    filtered_tokens.extend([word for word in tokens if word not in stopwords.words('english')])

                # Calculate word frequencies
                word_freq = Counter(filtered_tokens)

                # Sort the words by frequency in descending order
                sorted_word_freq = dict(sorted(word_freq.items(), key=lambda item: item[1], reverse=True))

                # Display only the top 100 frequencies
                number_of_words = top_of_word
                top_100_words = dict(list(sorted_word_freq.items())[:number_of_words])

                # Set Seaborn style
                sns.set(style="whitegrid")

                # Plot the Zipf distribution for the top 100 words using Seaborn
                plt.figure(figsize=(10, 6))
                ax = sns.barplot(x=list(top_100_words.values()), y=list(top_100_words.keys()), palette="viridis")
                ax.set(xlabel='Frequency', ylabel='Words')
                plt.title(f'Top {number_of_words} Zipf Distribution of Terms (Remove Stopwords)')
                plt.tight_layout()

                # Display the plot
                st.pyplot(plt)
            porter_algorithm  = st.sidebar.toggle("Porter’s algorithm", value=False)
            if porter_algorithm:

                # Remove stopwords and tokenize the text
                filtered_tokens = []
                for doc in documents:
                    tokens = clean_and_tokenize(doc)
                    filtered_tokens.extend([word for word in tokens if word not in stopwords.words('english')])

                # Apply Porter's stemming algorithm to the filtered tokens
                stemmer = PorterStemmer()
                stemmed_tokens = [stemmer.stem(word) for word in filtered_tokens]

                # Calculate word frequencies
                word_freq = Counter(stemmed_tokens)

                # Sort the words by frequency in descending order
                sorted_word_freq = dict(sorted(word_freq.items(), key=lambda item: item[1], reverse=True))

                # Display only the top 20 frequencies
                number_of_words = top_of_word
                top_words = dict(list(sorted_word_freq.items())[:number_of_words])

                # Set Seaborn style
                sns.set(style="whitegrid")

                # Plot the Zipf distribution for the top 20 stemmed words using Seaborn
                plt.figure(figsize=(10, 6))
                ax = sns.barplot(x=list(top_words.values()), y=list(top_words.keys()), palette="viridis")
                ax.set(xlabel='Frequency', ylabel='Words')
                plt.title(f'Top {number_of_words} Zipf Distribution of Terms (with Stopwords Removed and Porter Stemming)')
                plt.tight_layout()

                # Display the plot
                st.pyplot(plt)


                compare  = st.sidebar.toggle("compare the difference", value=False)
                if compare: 

                    # Remove stopwords and tokenize the text
                    filtered_tokens = []
                    for doc in documents:
                        tokens = clean_and_tokenize(doc)
                        filtered_tokens.extend([word for word in tokens if word not in stopwords.words('english')])

                    # Apply Porter's stemming algorithm to the filtered tokens
                    stemmer = PorterStemmer()
                    stemmed_tokens = [stemmer.stem(word) for word in filtered_tokens]

                    # Calculate word frequencies before stemming
                    word_freq_before = Counter(filtered_tokens)

                    # Calculate word frequencies after stemming
                    word_freq_after = Counter(stemmed_tokens)

                    # Sort the words by frequency in descending order
                    sorted_word_freq_before = dict(sorted(word_freq_before.items(), key=lambda item: item[1], reverse=True))
                    sorted_word_freq_after = dict(sorted(word_freq_after.items(), key=lambda item: item[1], reverse=True))

                    # Display only the top 20 frequencies
                    number_of_words = top_of_word
                    top_words_before = dict(list(sorted_word_freq_before.items())[:number_of_words])
                    top_words_after = dict(list(sorted_word_freq_after.items())[:number_of_words])

                    # Set Seaborn style
                    sns.set(style="whitegrid")

                    # Plot the Zipf distribution before and after Porter's stemming
                    plt.figure(figsize=(14, 6))
                    plt.subplot(1, 2, 1)
                    ax = sns.barplot(x=list(top_words_before.values()), y=list(top_words_before.keys()), palette="viridis")
                    ax.set(xlabel='Frequency', ylabel='Words')
                    plt.title(f'Top {number_of_words} Zipf Distribution of Terms (Before Stemming)')

                    plt.subplot(1, 2, 2)
                    ax = sns.barplot(x=list(top_words_after.values()), y=list(top_words_after.keys()), palette="viridis")
                    ax.set(xlabel='Frequency', ylabel='Words')
                    plt.title(f'Top {number_of_words} Zipf Distribution of Terms (After Stemming)')

                    plt.tight_layout()

                    # Display the plots
                    st.pyplot(plt)


            








