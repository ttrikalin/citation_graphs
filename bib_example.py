import os, re
import C, sel_isi , pubmed_functions, bib_to_graph


BASE_DIR = os.path.join(C.PUBMED_DATA_LOCATION, "citations") 

# this is the corpus title dictionary 
corpus = sel_isi.getPubMedPyData( C.PUBMED_DATA_LOCATION, C.DATA_LIMIT )

#lala = bib_to_graph.find_citedBy_in_corpus(BASE_DIR, corpus, "12598142")


D= bib_to_graph.make_citation_dict(BASE_DIR, corpus)

D2 = bib_to_graph.make_adjacency_matrix(D, "lala.txt", True)
print D2
