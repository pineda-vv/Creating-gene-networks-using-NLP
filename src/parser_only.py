from __future__ import division
from Bio import Medline
import pandas as pd
import string
import numpy as np


def abstract_dataframe(filename):
    """
    Input - file of Medline text output from Pubmed
    Output - DataFrame with PMID, Title for use in a psql db
    """
    pmid_ab_dict = medline_parser(filename)
    df = pd.DataFrame.from_dict(pmid_ab_dict, orient='index').reset_index()
    df.columns = ['pmid', 'title']
    df.to_csv('../data/pmid_titles_metabolism_5years.csv', index=False, index_label=False)



def medline_parser(filename):
    """extracts info from medline text file from pubmed"""
    pmid_abstract_dict = {}
    with open(filename) as handle:
        for record in Medline.parse(handle):
            if 'TI' in record.keys():
                pmid, title = record['PMID'], record['TI']
                pmid_abstract_dict[pmid] = title
        return pmid_abstract_dict

if __name__ == "__main__":
    """first file is PLOS Genetics abstract through June 2017"""
    # abstract_dataframe("../capstone_files/pubmed_result_medline.txt")
    # abstract_dataframe("../capstone_files/pubmed_result_plos_med.txt")
    # abstract_dataframe("../capstone_files/pubmed_result_plos_one.txt")
    # abstract_dataframe("../capstone_files/nature_genetics_all.txt")
    """genetics search term - filter reviews and last 5 years"""
    # abstract_dataframe("..capstone/genetics_search_reviews_5years.txt")
    # abstract_dataframe("../capstone_files/metabolism_5year_reviews.txt")
