import streamlit as st
from streamlit_option_menu import option_menu
import xml.etree.ElementTree as ET
import xml.etree.ElementTree as ET
from src.search_engine import parse_xml
from src.search_engine import search_and_highlight
import re

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
elif selected == "Frequency Analyst":
    st.image("image/analyst_title.png")

