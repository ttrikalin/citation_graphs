import os
from Bio import Entrez
# Tell Entrez who you are
Entrez.email = "ttrikalin@gmail.com"


def download_from_pubmed(search_term, filename, ret_number):
    """Downloads up to 2000 abstracts as an xml and saves it"""
    handle = Entrez.esearch(db="pubmed", term=search_term, retmode="xml", 
        RetMax=ret_number)
    
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
    authors = []
    for citation in citations:
        ids.append(citation["MedlineCitation"]["PMID"])
        names = []
        author_list = citation["MedlineCitation"]["Article"]["AuthorList"]
        for who in range(len(author_list)): 
            if author_list[who].keys()[0]=="LastName":
                lastname = author_list[who]["LastName"]
	        initials = author_list[who]["Initials"]
                new_name = lastname.lower() + "_" + initials.lower() 
	    else: 
	        newname = author_list[who]["CollectiveName"].lower()
	    names.append(new_name)
        authors.append(names)

    # close the handle!!! 
    in_handle.close()
    return (ids, authors)


