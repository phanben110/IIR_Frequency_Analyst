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
def download_pubmed(): 
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