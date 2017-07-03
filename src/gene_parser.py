import pandas as pd
import numpy as np
"""
input ==> dataframe with abstract and tokenized abstract
output ==> dataframe with gene pairs(or more) that can be set to
a dictionary for graph
"""
def generator():
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
        df_ucsc.columns = ['number', 'gene_name']
        gene_ucsc = set([str(name).lower() for name in df_ucsc["gene_name"] if len(str(name)) >1])
        return gene_ucsc

def gene_pair_finder(df, tokenized_col_name):
    """ input - dataframe, tokenized column name, set of gene names
        returns - converted dataframe with new column 'gene_pairs'
    """
    df['gene_pairs'] = df[tokenized_col_name].apply(if_gene_)
    return df

def if_gene_(df_column):
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

if __name__ == "__main__":
    generator()
