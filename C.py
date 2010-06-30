import C_login, os

PROJECT_NAME = 'Proton_beam'
DATA_LIMIT_H = 3200
DATA_LIMIT_L = 1600

if 10==10: 
    DEFAULT_SAVE_FILENAME = 'savedrecs.bib'
    SAVE_DIR ='C:\\Documents and Settings\\TTrikalinos\\Desktop\\test\\Downloads\\'
    DATA_DIR =  os.path.join('C:\\Documents and Settings\\TTrikalinos\\Desktop\\citation_graphs\\Results\\',PROJECT_NAME, 'citations')
    PUBMED_DATA_LOCATION = os.path.join('C:\\Documents and Settings\\TTrikalinos\\Desktop\\citation_graphs\\Results\\',PROJECT_NAME)
else:
    DEFAULT_SAVE_FILENAME = 'savedrecs.bib'
    SAVE_DIR = '/Users/tom/src/citation_graphs/temp' 
    DATA_DIR = os.path.join('/Users/tom/src/citation_graphs/Results/', PROJECT_NAME, 'citations')
    PUBMED_DATA_LOCATION = os.path.join('/Users/tom/src/citation_graphs/Results/', PROJECT_NAME)



# if the following two are identical and point to the login page I will login 
WEB_OF_KNOWLEDGE_URL_LOGIN = 'https://login.ezproxy.library.tufts.edu/login?auth=test&url=http://isiknowledge.com'
START_URL = WEB_OF_KNOWLEDGE_URL_LOGIN
IN_TUFTS_NETWORK = False 


# this is to set up a delay if we go too fast 
# Downloading approximately 140 per hour is not a big problem - when at home 
# but ~200 per hour at tufts makes it stop earlier...  so let's make silly delays... 
# 200 is 1 every 18 seconds
# 140 is 1 every 25 seconds (no delay)
# per paper set a delay of 7 sec
DELAY_IN_SECS = 7

# take a e.g 5 minute break every 50 papers
DELAY_BATCH_SIZE = 50
DELAY_PER_BATCH_IN_SECS = 5*60

# this is for tom 
LOGIN_USER = C_login.USER
LOGIN_PASS = C_login.PASS

