import os, re, pickle
import C, sel_isi , pubmed_functions, bib_to_graph, graph_functions
from igraph import * 

MY_XML = 'Pubmed.xml'
MY_GML = 'Graph.xml'
BASE_DIR = os.path.join(C.PUBMED_DATA_LOCATION, "citations") 

def main():
    # paths
    graph_path = os.path.join(C.PUBMED_DATA_LOCATION,'graph')
    data_path = os.path.join(C.PUBMED_DATA_LOCATION,'search')
    if (os.path.isdir(graph_path)==False ):
        os.makedirs(graph_path)
    if (os.path.isfile(os.path.join(data_path,MY_XML))==False ):
        print "There is no data xml!\n"
        exit ()

    # read the rcts and the observational studies
    # these wil be added as attributes to the Graph in a few lines below
    current_wd = os.getcwd()
    try:
        os.chdir(data_path)
        import Search
        reload(Search) # in case it has been edited and is being rerun in the same session 
        ALL_INDEX_PMIDS = []
        ALL_INDEX_PMIDS.extend(Search.rct_pmid_list)
        ALL_INDEX_PMIDS.extend(Search.obs_pmid_list)
        
    except:
        print 'Error in accessing "Search.py" in the directory\n'
        print os.path.join(data_path)
        exit ()
    os.chdir(current_wd)

    #try to read the graph if there is none - to avoind killing the process all the time
    if (os.path.isfile(os.path.join(graph_path,MY_GML))==False ):
    
        # this is the corpus title dictionary 
        corpus = sel_isi.getPubMedPyData( C.PUBMED_DATA_LOCATION)

        pmid_list, title_list, year_list = pubmed_functions.parse_pubmed_xml(os.path.join(data_path,MY_XML))
    
        D= bib_to_graph.make_citation_dict(BASE_DIR, corpus)
        node_dictionary, M = bib_to_graph.make_adjacency_matrix(D)
        G = Graph.Adjacency(M, mode="directed")
        G = graph_functions.add_standard_attributes_to_all_nodes(G, pmid_list, year_list, title_list, node_dictionary)
    
        G = graph_functions.set_attributes_for_specific_nodes(G, Search.rct_pmid_list, attribute_name='color', val_true='red', val_false='NA')
        G = graph_functions.set_attributes_for_specific_nodes(G, Search.obs_pmid_list, attribute_name='color', val_true='blue', val_false=None) # essential that val_false=None
    
        # add the rcts and obs list, the search string, date, size, and the topic title
        graph_functions.add_graph_attributes(G, ['rct_pmid_list', 'obs_pmid_list', 'search_string', 'search_date', 'corpus_size', 'topic_name'], [Search.rct_pmid_list, Search.obs_pmid_list, Search.search_string, Search.search_date, Search.corpus_size, Search.topic_name] )

        g = graph_functions.isolate_subcomponents_based_on_node_list(G, ALL_INDEX_PMIDS)

        if (os.path.isdir(graph_path)==False ):
            os.makedirs(graph_path)
            
        g.write(os.path.join(graph_path, MY_GML), format='graphml')
    
        #pickle objects to the Graph directory
        os.chdir(graph_path)
        P_data = '[P_data, D, node_dictionary, M, G, g]'
        P_list =[P_data, D, node_dictionary, M, G, g]
        #P_list = P_data
        f = open('Graph.pickle', 'wb')
        pickle.dump(P_list, f)
        f.close()
        print '[P_data, D, node_dictionary, M, G, g] pickled!'
        # return to the current directory
        os.chdir(current_wd)
    
        print 'Graph constructed succesfully; resides in \n', graph_path
        return(P_data, D, node_dictionary, M, G, g)

    else:
        #pickle objects to the Graph directory
        os.chdir(graph_path)
        f = open('Graph.pickle', 'rb')
        L= pickle.load(f)
        f.close()
        print '[P_data, D, node_dictionary, M, G, g] unpickled!'
        P_data = L[0]
        D = L[1]
        node_dictionary = L[2]
        M = L[3]
        G = L[4]
        g = L[5]
        # return to the current directory
        os.chdir(current_wd)
        return(P_data, D, node_dictionary, M, G, g)
        

P_data, D, node_dictionary, M, G, g = main()

#if (os.path.join(C.PUBMED_DATA_LOCATION,'graph', MY_GML)==False):
#    G,g = main()
#else:
#    P_data, D, node_dictionary, M, G, g = main()