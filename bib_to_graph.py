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
       
    print 'Running make_citation_dict...\n',
    
    D= {}
    i=0 # for the progress bar
    n= len(corpus_dict) # for the progress bar
    for index in corpus_dict:
        key, val = find_citedBy_in_corpus(base_dir, corpus_dict, index)
        D[key] = val
        # this is a progress report
        if (i%(n/10)==0):
            print '..',i/(n/10)*10,         
        i+=1

    print '\n',     
    return D


def write_adjacency_matrix_to_tabed_file(citation_dict, outfile, useheader=False):
    '''This prints out a text version of a citation matrix.
       Header
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
    if useheader ==True:
        header = 'row/col'
        empty_line = ''
        for col in range(n):
            header = header + '\t'+ D[col]
            empty_line += '\t0'
    
        f.write(header+ '\n')
    
    # write each row
    for row in range(n):
        if useheader == True: 
            myline = D[row]
        else:
            myline = None
        
        cited_by_list = citation_dict[D[row]]
        cited_by_string = "A" + "A".join(cited_by_list) + "A"
        
        if cited_by_string == '':   # if empty continue
            f.write(myline+empty_line+'\n')
            continue
        
        for col in range(n):   # not range(row, n) on purpose
                               # two papers could cite each other rarely!
            match = '0'
            if re.search("A" + D[col] + "A", cited_by_string):
                match ='1'
            myline = myline + "\t" + match
        
        f.write(myline + '\n')  # output line
    
    f.close()
    return D



def make_adjacency_matrix(citation_dict):
    '''This returns a list of matrix rows, which are lists themselves
       Also returns a Dictionary of the PMIDs that should be added 
       attributes to the vertices.
       ROW = INDEX
       COLUMN = CITED BY INDEX'''

    print 'Running make_adjacency matrix...\n',
    
    # this is an alphabetized dictionary of PMIDs
    pmid_list = citation_dict.keys()
    pmid_list.sort()
    n = len(pmid_list)
    D ={}
    for i, pmid in enumerate(pmid_list):
        D[i] = pmid

    M=[]   
    # write each row - now M is an upper triangular
    for row in range(n):
        row_list = []
        cited_by_list = citation_dict[D[row]]
        cited_by_string = "A"+"A".join(cited_by_list)+ "A"

        for col in range(n):   # not range(row, n) on purpose
                               # two papers could cite each other rarely!
            match = 0
            if re.search("A" + D[col] + "A", cited_by_string):
                match = 1
            row_list.append(match)
        
        M.append(row_list)
        # this is a progress report
        if (row%(n/10)==0):
            print '..',row/(n/10)*10, 
            
    print '\n',
    return D, M



