import os
from Bio import Entrez
# Tell Entrez who you are
Entrez.email = "ttrikalin@gmail.com"


def download_from_pubmed(search_term, filename, ret_number):
    """Downloads up to 2000 abstracts as an xml and saves it"""
    #handle = Entrez.esearch(db="pubmed", term=search_term, retmode="xml", 
    #    RetMax=ret_number)
    
    handle = Entrez.epost("pubmed", id=search_term)

    #record the id_list and then fetch the file into a local xml
    summary_record = Entrez.read(handle)
    id_list = summary_record["IdList"]

    if not os.path.isfile(filename):
        print "Downloading..."
        net_handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml" )
        out_handle = open(filename, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print "File \"" + filename + "\" saved!\n"
    elif os.path.isfile(filename):
        string= "File \"" + filename + "\" exists!\n"
        print string
    return



def parse_pubmed_xml(filename):
    """parse_pubmed_xml(filename): give a Pubmed xml and get back a list with the PMIDs 
       and a list of the lists of authors """
    in_handle = open(filename, "r")
    citations = Entrez.parse(in_handle)

    ids=[]
    titles = []
    for citation in citations:
        ids.append(citation["MedlineCitation"]["PMID"])
        titles.append(citation["MedlineCitation"]["Article"]["ArticleTitle"])

    # close the handle!!! 
    in_handle.close()
    return (ids, titles)




def prepare_search_files(ids_and_titles, rel_results_path, topic_name):

    '''takes a tuple of a list of pmids and titles and saves them in 
       the specified path, under "TI" '''
    current_path = os.path.abspath(".")
    topic_path  = os.path.join(current_path, rel_results_path, topic_name )
    
    if (os.path.isdir(topic_path)==False ):
        os.makedirs(topic_path)
    try:
        os.makedirs(topic_path+"/TI")
    except:
        print "Seems that this path exists:\n", topic_path+"/TI"

    os.chdir(topic_path+"/TI")
    
    for i, pmid in enumerate(ids_and_titles[0]):
        f = open(pmid, 'w')
        f.write(ids_and_titles[1][i].encode('ascii', 'replace'))
        f.close()

    os.chdir(current_path)

    print "Done! Data are in:\n", topic_path + "/TI"
