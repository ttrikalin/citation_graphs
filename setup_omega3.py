# Give a Pubmed XML and if needed edit the first line to automatically fix it :) 
# Then give the relative path to the Results folder in which the topic-specific 
# data will be saved, in a folder that has the PROJECT_NAME

# May be convenient to name the 

FILENAME = "omega3.xml"
REL_RES_PATH = "Results"
PROJECT_NAME = "omega3"

from pubmed_functions import * 

# get the xml if not already there 
# TO DO: automatically fix the filter bug
#download_from_pubmed(SEARCH, FILENAME, 5000)

# parse the xml into [pmids, titles, years]
pmids_titles_years = parse_pubmed_xml(FILENAME)

# organize saved files
prepare_search_files(pmids_titles_years, REL_RES_PATH, PROJECT_NAME)
