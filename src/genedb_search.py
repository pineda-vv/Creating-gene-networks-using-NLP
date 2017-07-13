import psycopg2 as pg
import networkx as nx
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from gephistreamer import graph
from gephistreamer import streamer
import numpy as np
import pandas as pd

def db_search():
    # Connect to an existing database
    conn = pg.connect("dbname=capstoneproject user=victorpineda")
    # Open a cursor to perform database operations
    cur = conn.cursor()
    gene = raw_input("Gene name :")
    query = "SELECT pmid FROM gene_pmid WHERE gene_name = '" + gene + "';"
    query2 = "SELECT gene_a, gene_b FROM combinations WHERE gene_a = '" + gene + "';"
    cur.execute(query2)
    ret = cur.fetchall()
    # gephi_streamer(ret)
    # ret_nx = {}
    # gene_list = []
    # for tup in ret:
    #     gene_list.append(tup[1])
    # ret_nx[gene] = gene_list
    feed = np.array(ret)
    df = pd.DataFrame(feed, columns=['gene_a', 'gene_b'])
    network_graph(df, 'gene_a', 'gene_b')

    # cur.fetchone()


    # Make the changes to the database persistent
    conn.commit()

    # Close communication with the database
    cur.close()
    conn.close()

# def gephi_streamer(tup_list):
#     streamer.GephiWS(hostname="localhost", port=8080, workspace="workspace0")
#     stream = streamer.Streamer(streamer.GephiWS())
#     for tup in tup_list:
#         stream.add_node(tup)
#     stream.commit()

def network_graph(df, col1, col2):
    """builds and visualizes networkx Graph for gene networks"""
    G = nx.from_pandas_dataframe(df, col1, col2)
    plt.figure(num=None, figsize=(16, 16), dpi=80, facecolor='w', edgecolor='c')
    nx.draw(G, with_labels=True, node_size=800, alpha=0.5, edge_color='c', cmap=plt.cm.Spectral)
    plt.show()

if __name__ == "__main__":
    db_search()

    topic_dict = {0: ' Plant Mechanisms',
             1: ' Clinical Trials',
             2: ' Cancer Therapeutics',
             3: ' Stem Cells and immune function',
             4: ' Energy Metabolism and Obesity',
             5: ' Protein-protein interactions',
             6: ' Receptor Activation',
             7: ' Diabetes and Obesity',
             8: ' Gene Expression Regulation',
             9: ' Stem Cell Growth/Survival',
             10: ' Transcriptional Mechanisms',
             11: ' Oxidative Stress',
             12: ' Brain/Blood Disorders',
             13: ' Bone Homeostasis',
             14: ' Methodology/Data Analysis',
             15: ' Inflammatory Response',
             16: ' Cardiovascular Biomarkers',
             17: ' Drug Targets',
             18: ' Membrane Lipids',
             19: ' Growth Signaling Pathway'}

    topic_number = {' Bone Homeostasis': 13,
             'Brain/Blood Disorders': 12,
             'Cancer Therapeutics': 2,
             'Cardiovascular Biomarkers': 16,
             'Clinical Trials': 1,
             'Diabetes and Obesity': 7,
             'Drug Targets': 17,
             'Energy Metabolism and Obesity': 4,
             'Gene Expression Regulation': 8,
             'Growth Signaling Pathway': 19,
             'Inflammatory Response': 15,
             'Membrane Lipids': 18,
             'Methodology/Data Analysis': 14,
             'Oxidative Stress': 11,
             'Plant Mechanisms': 0,
             'Protein-protein interactions': 5,
             'Receptor Activation': 6,
             'Stem Cell Growth/Survival': 9,
             'Stem Cells and immune function': 3,
             'Transcriptional Mechanisms': 10}
