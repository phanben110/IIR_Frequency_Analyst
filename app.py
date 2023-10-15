import streamlit as st
import xml.etree.ElementTree as ET
import re
import os 


st.set_page_config(page_title="Search PubMed Articles", page_icon="image/logo_csie2.png")
st.image("image/title_search.png")
st.sidebar.image("image/logo_NCKU.jpeg", use_column_width=True)


import xml.etree.ElementTree as ET

import xml.etree.ElementTree as ET

def parse_xml(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    data = []

    for article in root.findall('.//PubmedArticle'):
        # Initialize variables with empty strings
        pmid = ''
        title = ''
        abstract = ''
        journal_title = ''
        journal_issn = ''
        pubdate_year = ''
        pubdate_month = ''
        pubdate_day = ''
        author_list = []
        keyword_list = []

        pmid_element = article.find('.//PMID')
        if pmid_element is not None:
            pmid = pmid_element.text

        title_element = article.find('.//ArticleTitle')
        if title_element is not None:
            title = title_element.text

        abstract_element = article.find('.//Abstract/AbstractText')
        if abstract_element is not None:
            abstract = abstract_element.text

        # Extracting additional information
        journal_info = article.find('.//Journal')
        journal_title_element = journal_info.find('.//Title')
        if journal_title_element is not None:
            journal_title = journal_title_element.text

        journal_issn_element = journal_info.find('.//ISSN[@IssnType="Electronic"]')
        if journal_issn_element is not None:
            journal_issn = journal_issn_element.text

        pubdate = journal_info.find('.//PubDate')
        pubdate_year_element = pubdate.find('Year')
        pubdate_month_element = pubdate.find('Month')

        if pubdate_year_element is not None:
            pubdate_year = pubdate_year_element.text
        if pubdate_month_element is not None:
            pubdate_month = pubdate_month_element.text
        
        # Check if 'Day' element exists before accessing it
        pubdate_day_element = pubdate.find('Day')
        if pubdate_day_element is not None:
            pubdate_day = pubdate_day_element.text

        authors = article.find('.//AuthorList')
        if authors is not None:
            try:
                author_list = [f"{author.find('ForeName').text} {author.find('LastName').text}" for author in authors.findall('.//Author')]
            except AttributeError:
                author_list = []

        # Check if 'KeywordList' element exists before accessing it
        keyword_list_element = article.find('.//KeywordList[@Owner="NOTNLM"]')
        if keyword_list_element is not None:
            keyword_list = [keyword.text for keyword in keyword_list_element.findall('.//Keyword')]

        data.append({
            'PMID': pmid,
            'Title': title,
            'Journal Title': journal_title,
            'ISSN': journal_issn,
            'Publication Date': f"{pubdate_year}-{pubdate_month}-{pubdate_day}",
            'Abstract': abstract,
            'Authors': ', '.join(author_list),
            'Keywords': ', '.join(keyword_list)
        })

    return data



def search_and_highlight(article, search_term, case_sensitive=True):
    highlighted_fields = {}
    
    for key, value in article.items():
        flags = 0 if not case_sensitive else re.IGNORECASE
        try:
            highlighted_text = re.sub(
                fr'({re.escape(search_term)})',
                r'<span style="background-color: yellow">\1</span>',
                value,
                flags=flags,
            )
            if highlighted_text is not None:
                highlighted_fields[key] = highlighted_text
            else:
                highlighted_fields[key] = value
        except TypeError:
            # Handle the error by assigning the original value if 'value' is not a valid string
            highlighted_fields[key] = value

    return highlighted_fields


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
