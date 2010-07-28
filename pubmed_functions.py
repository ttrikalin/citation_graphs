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
       and a list of the lists of authors 
       all ordered by ascending alphanumerical pmid"""

    print 'Running parse_pubmed_xml...'
    
    in_handle = open(filename, "r")
    citations = Entrez.parse(in_handle)

    ids=[]
    titles = []
    years=[]
    
    for citation in citations:
        ids.append(citation["MedlineCitation"]["PMID"])
        titles.append(citation["MedlineCitation"]["Article"]["ArticleTitle"])
        #years.append(citation["PubMedData"]["History"]["PubMedPubDate"])
        list_of_year_dict = citation["PubmedData"]["History"]
        year = 3000 
        for d in list_of_year_dict:
            year = min(year, int(d["Year"]))
        years.append(year)
        
    # close the handle!!! 
    in_handle.close()

    D1={}
    D2={}
    for i in range(len(ids)):
        D1[i]=ids[i]
        D2[ids[i]]=i
    
    ids.sort()
    titles_sorted = []
    years_sorted = []
    for i in range(len(ids)):
        titles_sorted.append(titles[D2[ids[i]]] )
        years_sorted.append(years[D2[ids[i]]] )
    
    return (ids, titles_sorted, years_sorted)




def prepare_search_files(ids_titles_years, rel_results_path, topic_name):
    '''takes a tuple of a list of pmids, titles, years and saves them in 
       the specified path, under "TI" and "YR" respectively'''
    current_path = os.path.abspath(".")
    topic_path  = os.path.join(current_path, rel_results_path, topic_name )
    
    if (os.path.isdir(topic_path)==False ):
        os.makedirs(topic_path)
    try:
        os.makedirs(topic_path+"/TI")
    except:
        print "Seems that this path exists:\n", topic_path+"/TI"
    try:
        os.makedirs(topic_path+"/YR")
    except:
        print "Seems that this path exists:\n", topic_path+"/YR"

    #os.chdir(topic_path+"/TI")
    
    for i, pmid in enumerate(ids_titles_years[0]):
        f_ti = open(os.path.join(topic_path,'TI',pmid), 'w')
        f_ti.write(ids_titles_years[1][i].encode('ascii', 'replace'))
        f_ti.close()
        
        f_yr = open(os.path.join(topic_path,'YR',pmid), 'w')
        f_yr.write(str(ids_titles_years[2][i]))
        f_yr.close()

    #os.chdir(current_path)

    print "Done! Data are in:\n", topic_path + "[/TI, /YR]"
