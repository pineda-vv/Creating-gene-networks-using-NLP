from Bio import Medline
import pandas as pd
import string
import numpy as np
import nltk
import os
import glob
import unicodedata
from sklearn.feature_extraction.text import TfidfVectorizer, CountVectorizer
from sklearn.decomposition import NMF, LatentDirichletAllocation
from nltk.corpus import stopwords
from nltk.stem.snowball import SnowballStemmer
from nltk.stem import WordNetLemmatizer
from nltk.tokenize import word_tokenize
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import networkx as nx
import vtk

def abstract_dataframe(filename):
    """
    Input - file of Medline text output from Pubmed
    Output - DataFrame with PMID, Abstract, and Tokenized Abstract
    """
    pmid_ab_dict = medline_parser(filename)
    df = pd.DataFrame.from_dict(pmid_ab_dict, orient='index').reset_index()
    df.columns = ['PMID', 'Abstract']
    """tokenize abstract for NMF analysis"""
    df['tokenized_abs'] = df['Abstract'].apply(tokenize_abstract)
    df['gene_pairs'] = df['tokenized_abs'].apply(if_gene_)
    """create dictionary for networx_work"""
    gene_dict = {entry[0]:entry[1:] for entry in df['gene_pairs'] if entry != None}
    network_graph(gene_dict)
    topic_extraction(df, 'Abstract')


def medline_parser(filename):
    """input - medline text file from pubmed"""
    pmid_abstract_dict = {}
    with open(filename) as handle:
        for record in Medline.parse(handle):
            if 'AB' in record.keys():
                pmid, abstract = record['PMID'], record['AB']
                pmid_abstract_dict[pmid] = abstract
        return pmid_abstract_dict

def tokenize_abstract(abstract):
    """Abstract tokenizer"""
    stop = stopwords.words('english')
    clean_abstract = abstract.lower().translate(None, string.punctuation)
    tokens = word_tokenize(clean_abstract)
    stemmer = SnowballStemmer('english')
    stemmed_words = [stemmer.stem(word) for word in tokens]
    clean_words = [w for w in stemmed_words if w not in stopwords.words('english')]
    return ' '.join(clean_words)

def generator():
    """gene names generator"""
    filename_1 = '../capstone_files/ucsc_downloads/gene.txt'
    filename_2 = '../capstone_files/ucsc_downloads/geneSynonym.txt'
    gene_set_1 = gene_names(filename_1)
    gene_syn = gene_names(filename_2, complete=False)
    genes = gene_set_1 | gene_syn
    return genes

def gene_names(filepath, complete=True):
    """creates set from gene list file downloaded from UCSC Genome Browser"""
    if complete:
        df_ucsc = pd.read_csv(filepath, sep='\t', header=None)
        df_ucsc.columns = ['number', 'gene_name', 'locus_link', 'ref_seq_num', 'genbank', 'uniprot', 'taxon']
        gene_ucsc = set([str(name).lower() for name in df_ucsc["gene_name"] if len(str(name)) >1])
        return gene_ucsc
    else:
        df_syn = pd.read_csv(filepath, sep='\t', header=None)
        df_syn.columns = ['number', 'gene_name']
        gene_ucsc = set([str(name).lower() for name in df_syn["gene_name"] if len(str(name)) >1])
        return gene_ucsc

def if_gene_(abstract):
    genes = generator() # generates gene set
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
    if net_dict == None:
        net_dict = {}
    else:
        G = nx.from_dict_of_lists(net_dict)
    plt.figure(num=None, figsize=(100, 100), dpi=80, facecolor='w', edgecolor='k')
    nx.draw_networkx(G, node_shape='o', with_labels=True, alpha=0.5, edge_color='k', cmap=plt.cm.GnBu)
    plt.show()

def topic_extraction(df, col_name):
    tfidf_vectorizer = TfidfVectorizer(max_df=0.95, min_df=2,
                                       max_features=200,
                                       stop_words='english')
    tfidf = tfidf_vectorizer.fit_transform(df[col_name])

    tf_vectorizer = CountVectorizer(max_df=0.95, min_df=2,
                                    max_features=200,
                                    stop_words='english')
    tf = tf_vectorizer.fit_transform(df[col_name])
    nmf = NMF(n_components=20, random_state=1,
              alpha=.1, l1_ratio=.5).fit(tfidf)
    tfidf_feature_names = tfidf_vectorizer.get_feature_names()
    print_top_words(nmf, tfidf_feature_names)

    lda = LatentDirichletAllocation(n_topics=20, max_iter=5,
                                    learning_method='online',
                                    learning_offset=50.,
                                    random_state=0)
    lda.fit(tf)

    print("\nTopics in LDA model:")
    tf_feature_names = tf_vectorizer.get_feature_names()
    print_top_words(lda, tf_feature_names)

def print_top_words(model, feature_names, n_top_words=20):
    for topic_idx, topic in enumerate(model.components_):
        print("Topic #%d:" % topic_idx)
        print(" ".join([feature_names[i]
                        for i in topic.argsort()[:-n_top_words - 1:-1]]))
    print('End')


if __name__ == "__main__":
    # filename_list = ["../capstone_files/pubmed_result_medline.txt", "pubmed_result_plos_med.txt"]
    abstract_dataframe("../capstone_files/pubmed_result_medline.txt")
