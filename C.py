# the trik_test example!
import C_login


if 1==10: 
    DEFAULT_SAVE_FILENAME = 'savedrecs.bib'
    SAVE_DIR = 'C:\\Documents and Settings\\TTrikalinos\\Desktop\\test\\Downloads\\'
    DATA_DIR = 'C:\\Documents and Settings\\TTrikalinos\\Desktop\\citation_graphs\\Results\\omega3\\citations\\'
    PUBMED_DATA_LOCATION = 'C:\\Documents and Settings\\TTrikalinos\\Desktop\\citation_graphs\\Results\\omega3\\'
    DATA_LIMIT_H = None
    DATA_LIMIT_L = None
else:
    DEFAULT_SAVE_FILENAME = 'savedrecs.bib'
    SAVE_DIR = '/Users/tom/src/citation_graphs/temp' 
    DATA_DIR = '/Users/tom/src/citation_graphs/Results/omega3/citations' 
    PUBMED_DATA_LOCATION = '/Users/tom/src/citation_graphs/Results/omega3'
    DATA_LIMIT_H = None
    DATA_LIMIT_L = None



# if the following two are identical and point to the login page I will login 
WEB_OF_KNOWLEDGE_URL_LOGIN = 'https://login.ezproxy.library.tufts.edu/login?auth=test&url=http://isiknowledge.com'
START_URL = WEB_OF_KNOWLEDGE_URL_LOGIN
IN_TUFTS_NETWORK = False 

LOGIN_USER = C_login.USER
LOGIN_PASS = C_login.PASS

