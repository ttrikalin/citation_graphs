import os, re
import C, sel_isi , pubmed_functions, bib_to_graph
from igraph import * 


BASE_DIR = os.path.join(C.PUBMED_DATA_LOCATION, "citations") 

# this is the corpus title dictionary 
corpus = sel_isi.getPubMedPyData( C.PUBMED_DATA_LOCATION)

#lala = bib_to_graph.find_citedBy_in_corpus(BASE_DIR, corpus, "12598142")

years = pubmed_functions.parse_pubmed_xml('Trikalinos.xml')[2]

D= bib_to_graph.make_citation_dict(BASE_DIR, corpus)

pmid, M = bib_to_graph.make_adjacency_matrix(D)

G = Graph.Adjacency(M, mode="directed")



for i, v in enumerate(G.vs):
    v["pmid"] = pmid[i]
    v['year'] = years[i]
    v['design'] = None
    v['inma'] = True

print G.vs[3]
G.es