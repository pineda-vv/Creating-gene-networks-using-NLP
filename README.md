# GAIN - Gene Association/Interaction Network
## Creating a gene-association network from text taken from abstracts of research articles in Pubmed
### A capstone project for the Galvanize Data Science Immersive

Gene Network Graph derived from abstracts taken from five years of metabolism review articles
![Alt text](https://github.com/pineda-vv/Creating-gene-networks-using-NLP/blob/master/data/metab_with_labels.png)

#### Business Understanding
---
An organism’s set of genes provides the basic blueprint for the characteristic features that define what that organism.  The interaction between the products of these genes dictates the cellular and organismal response to stimuli, environmetal changes, and from input from other cells or organisms. As an example, a complex network of the conserved genes dictates the process of embryonic development in different multi-cellular organisms.   Our understanding of these networks emerged from a multitude of experimental evidence and each layer is added constantly as our knowledge base increases.

The question I asked is the following.  Can I build a gene network without prior knowledge of biology?  Can I mine the text of scientific abstracts to create such a network?  And, what is the utility of this endeavor?  For the last question, one can imagine a relative science novice looking at the created gene network as a starting point to investigating a gene and its functional partners.  Another application may be for a company that has an interest with a specific disease-causing gene.  Finding the gene-pairs/interactors for the gene of interest may point to a new mechanistic pathway that may be an important pharmacologic target for the disease

#### Data Understanding
---
###### NLM/Pubmed

There are literally millions of published articles annotated in Pubmed which are readily available.  To not be overwhelmed with too much information, I started my analysis with an abstract list returned from a query for the keyword "metabolism" and filtering for just review articles and only from the last 5 years.  This analysis provides the proof of concept and results in a deliverable minimum viable product in the time that is allotted for the capstone project.  I would then expand the project to include abstracts from several years, which would then include most of the papers in genetics and biology after the human genome was mostly annotated.   Since some of the journals that are accessible by Pubmed cater to a wide audience (Science, Nature, Proceedings of the National Academy of Sciences), I will focus on just review articles during my data collection and/or subsequent analysis.

#### Data Preparation
---
The data pipeline for the text analysis begins with a Pubmed download of a formatted text file (Medline format).   The python module Biopython allows for parsing of several components of this text, treating it as a dictionary file.  I extracted the Pubmed ID (PMID), the title (TI), and the body of the abstract (AB).  The list of genes used to query the abstract text is built from two online repositories - the UCSC Genome Browser and Uniprot.  

#### Model
---
The abstracts extracted from pubmed were analyzed in two ways.  Firstly, each abstract went through a tokenizer (NLTK) and these tokenized abstracts were used to mine for the gene pairs.  Secondly, the abstracts were clustered into twenty categories mined from topic analysis using non-negative matrix factorization (NMF).  Twenty categories were chosen as an optimal number of categories/topics that had very little if any overlap.  

#### Deployment
---
The gene pairing and NMF analysis results were pooled and stored in a Postgres database that is used by a search engine that can return a networkx Graph for a gene of interest, titles of associated articles, and as a metric, a plot of the topic distribution for the articles.  After the dataset is expanded, a web-based search engine will be deployed that will take a gene name in any form and return a graph and article list.

Example graph centered on the brca2 gene
![Alt text](https://github.com/pineda-vv/Creating-gene-networks-using-NLP/blob/master/data/brca2_new.png)
