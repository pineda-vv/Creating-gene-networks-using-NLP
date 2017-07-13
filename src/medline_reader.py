from __future__ import division
from Bio import Medline
import pandas as pd
import string
import numpy as np
import multiprocessing as mp
import nltk
import os
import glob
import unicodedata
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import TfidfVectorizer, CountVectorizer
from sklearn.decomposition import NMF, LatentDirichletAllocation
from nltk.corpus import stopwords
from nltk.stem.snowball import SnowballStemmer
from nltk.stem import WordNetLemmatizer
from nltk.tokenize import word_tokenize
import matplotlib.pyplot as plt
import networkx as nx

def abstract_dataframe(filename):
    """
    Input - file of Medline text output from Pubmed
    Output - DataFrame with PMID, Abstract, and Tokenized Abstract
    """
    pmid_ab_dict = medline_parser(filename)
    df = pd.DataFrame.from_dict(pmid_ab_dict, orient='index').reset_index()
    df.columns = ['PMID', 'Abstract']
    """
    Parallelized tokenizer and gene pairs functions gene-network analysis.
    Original code used are commented out
    """
    df = parallel_tokenizer(df)
    df = parallel_genepairs(df)
    """create dictionary for networx_work"""
    df = topic_extraction(df, 'Abstract') # after topic extraction
    df.to_csv('metabolism_5years_tokenized.csv')
    gene_dict = {entry[0]:entry[1:] for entry in df['gene_pairs'] if entry != None}
    network_graph(gene_dict)


def medline_parser(filename):
    """extracts info from medline text file from pubmed"""
    pmid_abstract_dict = {}
    with open(filename) as handle:
        for record in Medline.parse(handle):
            if 'AB' in record.keys():
                pmid, abstract = record['PMID'], record['AB']
                pmid_abstract_dict[pmid] = abstract
        return pmid_abstract_dict

def parallel_tokenizer(df):
    """parallelizes the tokenizer function"""
    pool = mp.Pool(processes=4)
    df['tokenized_abs'] = pool.map(_tokenize_abstract, df['Abstract'])
    pool.terminate()
    return df

def parallel_genepairs(df):
    """parallelizes the gene pairing"""
    pool = mp.Pool(processes=4)
    df['gene_pairs'] = pool.map(_if_gene, df['tokenized_abs'])
    pool.terminate()
    return df

def _tokenize_abstract(abstract):
    """Abstract tokenizer"""
    stop = stopwords.words('english')
    clean_abstract = abstract.lower().translate(None, string.punctuation)
    tokens = word_tokenize(clean_abstract)
    stemmer = SnowballStemmer('english')
    stemmed_words = [stemmer.stem(word) for word in tokens]
    clean_words = [w for w in stemmed_words if w not in stopwords.words('english')]
    return ' '.join(clean_words)

def _generator():
    """gene names generator from two files"""
    filename_1 = 'gene.txt'
    filename_2 = 'geneSynonym.txt'
    gene_set_1 = gene_names(filename_1)
    gene_syn = gene_names(filename_2, complete=False)
    genes = gene_set_1 | gene_syn
    return genes

def gene_names(filepath, complete=True):
    """creates a set from 2 gene list files downloaded from the UCSC Genome Browser"""
    if complete:
        df_ucsc = pd.read_csv(filepath, sep='\t', header=None)
        df_ucsc.columns = (
                ['number', 'gene_name', 'locus_link',
                 'ref_seq_num', 'genbank', 'uniprot', 'taxon']
            )
        gene_ucsc = set(
                [str(name).lower() for name in df_ucsc["gene_name"]
                if len(str(name)) >1]
            )
        return gene_ucsc
    else:
        df_syn = pd.read_csv(filepath, sep='\t', header=None)
        df_syn.columns = ['number', 'gene_name']
        gene_ucsc = set(
                [str(name).lower() for name in df_syn["gene_name"]
                if len(str(name)) >1]
            )
        return gene_ucsc

def _if_gene(abstract):
    genes = _generator() # generates gene set
    gene_pairs = set()
    for item in abstract.split():
        word = item.strip('.').strip(',').lower()
        if word in genes:
            gene_pairs.add(word)
    if len(gene_pairs) > 1:
        return list(gene_pairs)
    else:
        return None

def network_graph(net_dict=None):
    """builds and visualizes networkx Graph for gene networks"""
    if net_dict == None:
        net_dict = {}
    else:
        G = nx.from_dict_of_lists(net_dict)
    plt.figure(num=None, figsize=(30, 30), dpi=80, facecolor='w', edgecolor='c')
    nx.draw_networkx(G, with_labels=True, alpha=0.5, edge_color='c', cmap=plt.cm.GnBu)
    plt.savefig("metabolism_5years.png", bbox_inches='tight')

def topic_extraction(df, col_name):
    """
    Two algorithms for topic extraction -
    NMF and LatentDirichletAllocation(LDA).
    Need to tie in with K-means clustering
    """
    tfidf_vectorizer = TfidfVectorizer(max_df=0.95, min_df=2,
                                       max_features=200,
                                       stop_words='english')
    tfidf = tfidf_vectorizer.fit_transform(df[col_name])

    tf_vectorizer = CountVectorizer(max_df=0.95, min_df=2,
                                    max_features=200,
                                    stop_words='english')
    tf = tf_vectorizer.fit_transform(df[col_name])
    nmf = NMF(n_components=20, random_state=1,
            alpha=.1, l1_ratio=.5)
    tfidf_feature_names = tfidf_vectorizer.get_feature_names()
    nmf_w = nmf.fit_transform(tfidf)
    nmf_h = nmf.components_
    df['labels'] = nmf_w.argmax(axis=1) # this was the right code to get labels/clusters


    print("\nTopics in NMF model:")
    print_top_words(nmf, tfidf_feature_names)


    lda = LatentDirichletAllocation(n_topics=20, max_iter=5,
                                learning_method='online',
                                learning_offset=50.,
                                random_state=0,
                                n_jobs=-1)
    lda.fit(tf)
    df['perplexity'] = lda.perplexity(tf)
    print("\nTopics in LDA model:")
    tf_feature_names = tf_vectorizer.get_feature_names()
    print_top_words(lda, tf_feature_names)
    return df

def print_top_words(model, feature_names, n_top_words=20):
    for topic_idx, topic in enumerate(model.components_):
        print("Topic #%d:" % topic_idx)
        print(" ".join([feature_names[i]
                        for i in topic.argsort()[:-n_top_words - 1:-1]]))
    print('End')


if __name__ == "__main__":
    """first file is PLOS Genetics abstract through June 2017"""
    # abstract_dataframe("../capstone_files/pubmed_result_medline.txt")
    # abstract_dataframe("../capstone_files/pubmed_result_plos_med.txt")
    # abstract_dataframe("../capstone_files/pubmed_result_plos_one.txt")
    # abstract_dataframe("../capstone_files/nature_genetics_all.txt")
    """genetics search term - filter reviews and last 5 years"""
    # abstract_dataframe("genetics_search_reviews_5years.txt")
    abstract_dataframe("metabolism_5year_reviews.txt")
