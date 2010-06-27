import os, re 
import C

def get_titles_from_bib(base_dir, pmid):
    '''Assumes that the name of the bib will be the pmid
       and returns the pmid and the list of titles.
       base_dir should be: 
       os.path.join(C.PUBMED_DATA_LOCATION, "citations") '''

    bibfile = os.path.join(base_dir, pmid+".bib")
    f = open(bibfile, "r")
    data = f.readlines()
    f.close()

    # join to form a single string 
    data2 = "".join(data)

    # remove newline characters and long whitespace
    data2 = re.sub("\n", " ", data2).strip()
    data2 = re.sub("[ ]{2,}", " ", data2).strip()

    # match the titles 
    titles = re.findall("Title = {{(.*?)}" , data2) 

    return pmid, titles 


def find_citedBy_in_corpus(base_dir, corpus_dict, index_pmid):
    '''Takes a dictionary for the corpus, the index pmid and the location
       of the bib_files and returns the pmids that are cited by the index pmid
       base_dir should be: 
       os.path.join(C.PUBMED_DATA_LOCATION, "citations") '''

    candidates = get_titles_from_bib(base_dir, index_pmid)[1]
    
    cited_pmid = []
    if len(candidates) == 0:
        return index_pmid, cited_pmid

    #iterate over corpus dict
    for ti in corpus_dict: 

        #iterate over cited titles
        for cited_ti in candidates:
    
            #check for match -- enforse ascii and remove brackets of non-english titles
            w1 = corpus_dict[ti]['title'].encode('ascii', 'replace').strip().lower().rstrip('].').lstrip('[')
            w2 = cited_ti.encode('ascii','replace').strip().lower().rstrip('].').lstrip('[') 

            if (w1==w2) & (w1!=""):
                cited_pmid.append(ti)

    #done
    return index_pmid, cited_pmid



def make_citation_dict(base_dir, corpus_dict):
    '''Takes a dictionary for the corpus, the index pmid and the location
       of the bib_files and returns the pmids that are cited by the index pmid
       base_dir should be: 
       os.path.join(C.PUBMED_DATA_LOCATION, "citations") '''

    D= {}
    for index in corpus_dict: 
        key, val = find_citedBy_in_corpus(base_dir, corpus_dict, index)
        D[key] = val

    return D


def make_adjacency_matrix_to_file(citation_dict, outfile, brief):
    '''This prints out a text version of a citation matrix that is either 
       brief or complete. 
       ROW = INDEX
       COLUMN = CITED BY INDEX'''


    
    # this is alphabetized dictionary of PMIDs
    list = citation_dict.keys()
    list.sort()
    n = len(list)
    D ={}
    for i, pmid in enumerate(list):
        D[i] = pmid
    

    # the file handle 
    f= open(outfile, 'w')
    
    # write the header --
    # and create and empty line 
    header = 'row/col' 
    empty_line = ''
    for col in range(n):
        header = header + '\t'+ D[col]
        empty_line += '\t0'

    f.write(header+ '\n')
        
    # write each row 
    for row in range(n):
        myline = D[row]

        cited_by_list = citation_dict[D[row]]
        cited_by_string = "|".join(cited_by_list)

        if cited_by_string == '':   # if empty continue 
            f.write(myline+empty_line+'\n')
            continue 

        for col in range(n):   # not range(row, n) on purpose
                               # two papers could cite each other rarely!
            match = '0' 
            if re.search(D[col], cited_by_string):
                match ='1'
            myline = myline + "\t" + match 

        f.write(myline + '\n')  # output line 
    
    f.close()
    return D
