#SEARCH = '(("dietary fats, unsaturated" [MH] OR "fish oils" [MH] OR "omega 3" OR "linolenic" OR "docosahexaenoic" OR "eicosapentaenoic" OR "Fish oil" OR "Fatty acids" OR "Fish consumption" OR "omega3") AND (cardiovascular disease OR mortality)) AND (cohort OR prospective trial OR clinical trial OR controlled trial OR longitudinal study OR prospective study OR randomised OR randomized OR Prevenzione)'
SEARCH = '(("dietary fats, unsaturated" [MeSH Terms] OR "fish oils" [MeSH Terms] OR "omega 3" [Text Word] OR "linolenic" [Text Word] OR "docosahexaenoic" [Text Word] OR "eicosapentaenoic" [Text Word] OR "Fish oil" [Text Word] OR "Fatty acids" [Text Word] OR "Fish consumption" [Text Word] OR "omega3" [Text Word]) AND (cardiovascular disease OR mortality)) AND (cohort OR prospective trial OR clinical trial OR controlled trial OR longitudinal study OR prospective study OR randomised OR randomized OR Prevenzione)'
FILENAME = "omega3.xml"
REL_RES_PATH = "Results"
TOPIC_NAME = "omega3"

from pubmed_functions import * 

# get the xml if not already there 
# TO DO: automatically fix the filter bug
#download_from_pubmed(SEARCH, FILENAME, 5000)

# parse the xml 
data = parse_pubmed_xml(FILENAME)

# organize saved files
prepare_search_files(data, REL_RES_PATH, TOPIC_NAME)
