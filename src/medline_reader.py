from Bio import Medline
import pandas as pd
import string
import numpy as np
import nltk
import os
import glob
import unicodedata
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.feature_extraction.text import CountVectorizer
from nltk.corpus import stopwords
from nltk.stem.snowball import SnowballStemmer
from nltk.stem import WordNetLemmatizer
import string
from nltk.tokenize import word_tokenize
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity

def abstract_dataframe(filename, output_name):
    pmid_ab_dict = medline_parser(filename)
    df = pd.DataFrame.from_dict(pmid_ab_dict, orient='index').reset_index()
    df.columns = ['PMID', 'Abstract']
    df['tokenized_abs'] = df['Abstract'].apply(tokenize_abstract)
    df.to_csv(output_name, sep='\t')

def medline_parser(filename):
    """input - downloaded medline text file from pubmed"""
    pmid_abstract_dict = {}
    with open(filename) as handle:
        for record in Medline.parse(handle):
            if 'AB' in record.keys():
                pmid, abstract = record['PMID'], record['AB']
                pmid_abstract_dict[pmid] = abstract
        return pmid_abstract_dict

def dictionary_compiler(dictionary):
    tokenized_dict = {}
    for k, v in abstract_dict.iteritems():
        token = tokenize_abstract(v)
        tokenized_dict[k] = token


def tokenize_abstract(abstract):
    stop = stopwords.words('english')
    clean_abstract = abstract.lower().translate(None, string.punctuation)
    tokens = word_tokenize(clean_abstract)
    stemmer = SnowballStemmer('english')
    stemmed_words = [stemmer.stem(word) for word in tokens]
    clean_words = [w for w in stemmed_words if w not in stopwords.words('english')]
    return ' '.join(clean_words)


if __name__ == "__main__":
    filename_list = ["../capstone_files/pubmed_result_medline.txt", "pubmed_result_plos_med.txt"]
    abstract_dataframe(filename_list[0],'plos_genetics.csv')
