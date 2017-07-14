import psycopg2 as pg
import networkx as nx
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import seaborn as sns
from gephistreamer import graph
from gephistreamer import streamer
import numpy as np
import pandas as pd
import pygraphviz
from networkx.drawing.nx_agraph import graphviz_layout

def db_search():
    # Connect to an existing database
    conn = pg.connect("dbname=capstoneproject user=victorpineda")
    # Open a cursor to perform database operations
    cur = conn.cursor()
    query0 = "SELECT gene_name from gene_pmid;"
    cur.execute(query0)
    gene_list = cur.fetchall()
    gene_set = set(tup[0] for tup in gene_list)

    gene_quest = False
    while gene_quest == False:
        gene_input = raw_input("Gene name :")
        gene = gene_input.lower()
        if gene in gene_set:
            gene_quest = True
    query = "SELECT pmid FROM gene_pmid WHERE gene_name = '" + gene + "';"
    cur.execute(query)
    pmid_list = cur.fetchall()
    titles = []
    for pmid in pmid_list:
        art_query = "SELECT title from articles where pmid = {}".format(pmid[0])
        cur.execute(art_query)
        title = cur.fetchone()
        pm_title = (int(pmid[0]), title[0])
        titles.append(pm_title)
    print "Five relevant reviews"
    print "=" * 20
    if len(titles) < 5:
        for tup in titles:
            print tup
    else:
        for tup in titles[:5]:
            print tup
    query2 = "SELECT gene_a, gene_b FROM combinations WHERE gene_a = '" + gene + "';"
    cur.execute(query2)
    ret = cur.fetchall()
    ret_dict = Counter(ret)
    ret_arr = [(k[0], k[1], v) for k, v in ret_dict.items()]
    feed = np.array(ret_arr)
    df = pd.DataFrame(feed, columns=['gene_a', 'gene_b', 'count'])
    # network_graph(df, 'gene_a', 'gene_b', 'count')
    weighted_network_graph(ret_dict)
    # cur.fetchone()


    # Close communication with the database
    cur.close()
    conn.close()

def network_graph(df, col1, col2, col3):
    """builds and visualizes networkx Graph for gene networks"""
    G = nx.from_pandas_dataframe(df, col1, col2, edge_attr=col3)
    plt.figure(num=None, figsize=(16, 16), dpi=80, facecolor='w', edgecolor='c')
    nx.draw(G, pos=graphviz_layout(G), with_labels=True, node_size=1600, cmap=plt.cm.Spectral,
        node_color=range(len(G)),
        prog='dot', font_color='k', font_weight='bold')
    nx_draw_edges(G)
    plt.show()

def weighted_network_graph(counter_dict):
    G = nx.Graph()
    for k, v in counter_dict.items():
        G.add_edge(k[0], k[1], weight = v)
    elarge = [(u,v) for (u,v,d) in G.edges(data=True) if d['weight'] > 5]
    emedium = [(u,v) for (u,v,d) in G.edges(data=True) if d['weight'] > 1 and d['weight'] < 4]
    esmall = [(u,v) for (u,v,d) in G.edges(data=True) if d['weight'] == 1]
    pos = graphviz_layout(G)
    plt.figure(num=None, figsize=(16, 16), dpi=80, facecolor='w', edgecolor='c')
    nx.draw_networkx_nodes(G, pos
                    , node_size=1600, with_labels=True
                    , node_color=range(len(G))
                    , cmap=plt.cm.Spectral
                )
    nx.draw_networkx_edges(G, pos
                    , edgelist=elarge
                    , width=5
                )
    nx.draw_networkx_edges(G, pos
                    , edgelist=emedium
                    , width=2
                    , edge_color='b'
                    , alpha=0.7
                )
    nx.draw_networkx_edges(G, pos
                    , edgelist=esmall
                    , width=1,alpha=0.5
                    , edge_color='r'
                    , style='dashed'
                )
    nx.draw_networkx_labels(G,pos, font_size=12
                    , font_family='sans-serif'
                    , font_color='k'
                    , font_weight='bold'
                )
    plt.axis('off')
    # plt.savefig("../data/weighted_ppard.png") # save as png
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
