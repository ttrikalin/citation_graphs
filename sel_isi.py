import cookielib
from collections import defaultdict
import os, sys
from sys import exit, stdin, stdout, stderr
from BeautifulSoup import BeautifulSoup
from selenium import selenium
import time
import C

class Browser:

	def __init__(self):
		self.PAGE_LOAD_TIMEOUT = 1000
		self.url = 'http://apps.isiknowledge.com'
		self.br = selenium("localhost", 4444, "*chrome", "http://apps.isiknowledge.com")
		#self.br = selenium("localhost", 4444, "*googlechrome /Applications/Google Chrome.app/Contents/MacOS/Google Chrome", "http://apps.isiknowledge.com")
		self.FILE_SAVING_TIME = 1.0   # seconds to wait before attempting to read downloaded files
		self.SKIP_TIME = 0.5    # seconds to wait before re-attempting to read downloaded files

        def create_stub_file(self, pmid):
            f = open(  os.sep.join( [C.DATA_DIR, pmid+'.bib'] ) ,  'w'  )
            f.close()

        def stop(self):
            self.br.stop()

	def run_article_TEST(self, title, pmid):
                self.robust_start()
                title = title.lstrip('['); title = title.rstrip(']');
		cb_link, numCiting = self.get_citedByLink(title)
                if cb_link is None:
                    self.br.stop(); return;
		print_debug( cb_link + ' with ' + str(numCiting) + ' citations\n'  )
		self.get_citedBy_data(cb_link, numCiting, pmid)
		self.br.stop()

	def run_article(self, article):

                self.robust_start()                                             # Start Browser session for this article (required)
            
		try:
			pmid = article['pmid']
			title = article['title']
		except:
                        self.create_stub_file(pmid)                             # Article doesn't have a title. Creating an empty file is OK in this case.
                        self.br.stop()
			return

		titleSearchString = title.rstrip(']').lstrip('[')

		cb_link, numCiting = self.get_citedByLink( titleSearchString )            # Get the cited by url (link)
		
                if cb_link is None:
                    self.create_stub_file(pmid)         # No citedBy link means no citations, so create an empty citations file.
                    self.br.stop()
                    return
                
		self.get_citedBy_data(cb_link, numCiting, pmid)

		self.br.stop()                                             # Stop Browser session for this article (required)


	def get_citedByLink(self, titleSearchString ):

		GOLD_Title = titleSearchString
		
		########## TESTING ONLY ############
                if in_test(): titleSearchString = 'AMP-activated protein kinase phosphorylation of endothelial NO synthase.'
                ####################################

                succ = self.robust_openAndWait(self.url)
                if not succ: return None, None
			
		self.br.type( 'value(input1)' ,  titleSearchString )
		self.br.select( 'value(select1)' , 'Title' )	# Choose 'Title' from the dropdown selecting the string search target field
		self.br.click( """//input[@title='Search']""")	# 

                wait_succ = self.robust_wait()
                if not wait_succ: return None, None

		html_u = self.br.get_html_source()
		html = html_u.encode('ascii', 'replace')
		soup = BeautifulSoup(html)
		results = soup.findAll('a', {'class':'smallV110'})
		found = False
		for r in results:
			this_title =  str(r.contents[0])
                        print_debug( "[ Do Titles Match? ] : " + this_title )
			if this_title.strip().lower().rstrip('].').lstrip('[') == GOLD_Title.strip().lower().rstrip('].').lstrip('[') :
				print_debug( ' [ YES ... FOUND A MATCH ]\n' )
				try:
                                        citedLink = r.findNext('a',title='View all of the articles that cite this one')['href']
                                        numCiting = r.findNext('a',title='View all of the articles that cite this one').contents
                                        numCiting = int(str(numCiting[0]).replace(',',''))
                                        found = True
                                        break
                                except:
                                       continue
                        else:
                               print_debug( ' [ no ]\n' ) 
		if found: return citedLink, numCiting
		else:     return None, None


	def get_citedBy_data(self, cb_link, numCiting, pmid):

            ow_succ = self.robust_openAndWait( self.url + cb_link )
            if not ow_succ:
                ### self.create_stub_file(pmid)                 ### Web error occurred, dont assume article has no citations, dont create empty file
                return
                
            self.br.select('pageSize', "Show 50 per page" )     ### Max is 50 per page, url tricks wont increase this either
            
            w_succ = self.robust_wait()
            if not w_succ:
                #self.create_stub_file(pmid)                    ### Web error occurred, dont assume article has no citations, dont create empty file
                return
                
            self.do_save( pmid , 1 )        ## Save this 1st page
            
            # Build new page link prefix, and find the number of pages to be processed
            ar = cb_link.split('&'); newLinkPrefix = '&'.join( [ '/summary.do?product=UA', 'qid=3', ar[1], ar[2] ] );
            numPages = self.calc_num_pages(numCiting)

            # Save pages 2 -> Last
            for pageNum in xrange(2, numPages+1):
                newLink = self.url + newLinkPrefix + '&formValue(summary_mode)=CitingArticles' + '&page=' + str(pageNum)
                self.br.stop()      # Get around XSS security provided by the browser
                time.sleep(0.2)
                self.robust_start()      # Get around XSS security provided by the browser
                succ = self.robust_openAndWait( newLink )
                if not succ: continue
                self.do_save(pmid, pageNum)

	def robust_start(self):
            for i in xrange(300):
                try:
                    self.br.start()
                    return True
                except:
                    time.sleep(0.2)
            return False

	def robust_open(self,url):
            for i in xrange(300):
                try:
                    self.br.open(url)
                    return True
                except:
                    time.sleep(0.2)
            return False

        def robust_wait(self):
            for i in xrange(300):
                try:
                    self.br.wait_for_page_to_load( self.PAGE_LOAD_TIMEOUT )
                    return True
                except:
                    time.sleep(0.2)
            return False
        

        def robust_openAndWait(self, url):

            open_succ = self.robust_open(url)
            if not open_succ: return False

            wait_succ = self.robust_wait()
            if not wait_succ: return False

            return True
                
        def file_append( self, src, out, pageNum ):

            r1 = True
            r2 = True
            read = False

            # Make sure this particular file ('ABCD.txt') is not there. We need this reserved for a rename-flip file access test.
            try:
                os.remove( C.SAVE_DIR+'ABCD.txt' )
            except:
                pass
                      
            for i in xrange(int(30/self.SKIP_TIME)):     # Keep trying for 30 seconds

                # Sections r1 and r2 do a rename-flip to see if the OS will allow us to mangle this file yet.
                # If the OS lets us mangle it, then we assume that the file has downloaded fully and it's
                # now OK to read its contents
                if r1:
                    try:
                        os.rename( src , C.SAVE_DIR+'ABCD.txt' )
                        r1 = False
                    except:
                        time.sleep( self.SKIP_TIME )
                        continue
                
                if r2:
                    try:
                        os.rename( C.SAVE_DIR+'ABCD.txt', src )
                        r2 = False
                    except:
                        time.sleep( self.SKIP_TIME )
                        continue

                try:
                    dummy = open( src , 'a' )   # Draw an exception if we're not allowed to write to this file yet
                    dummy.close()
                except:
                    time.sleep( self.SKIP_TIME )

                try:
                    src_fi = open(src, 'r')             # read the file's contents
                    data = src_fi.readlines()
                    src_fi.close()
                    read = True
                    break
                except:
                    time.sleep( self.SKIP_TIME )        # Wait SKIP_TIME seconds
                    
                
            if not read: return False

            if pageNum == 1:    WRITE_MODE = 'w'
            else:               WRITE_MODE = 'a'

            if not os.path.isdir( C.DATA_DIR ):
                os.mkdir( C.DATA_DIR )
                
            out_f = open(out, WRITE_MODE)
            for l in data:
                out_f.write(l)
            out_f.close()
            return True


        def kill_recs_file(self):
            try:    os.remove(C.SAVE_DIR + C.DEFAULT_SAVE_FILENAME)
            except: pass

        def do_append(self, pmid, pageNum):
            appended = False
            # Try 25 times to append. Each append may take 30sec.
            for i in xrange(25):
                appended = self.file_append( C.SAVE_DIR+C.DEFAULT_SAVE_FILENAME, C.DATA_DIR+pmid+'.bib', pageNum )
                if appended: break
            self.kill_recs_file()

        def do_save(self, pmid, pageNum):
            self.kill_recs_file()
            self.setup_save_options()
            self.br.click("save")
            # Check to see if the newly created file is there yet (old one gets deleted every time, hence we can use this as a valid test)
            while not os.path.isfile( C.SAVE_DIR+C.DEFAULT_SAVE_FILENAME ):
                time.sleep(1.0)
            self.do_append(pmid, pageNum)

        def go_to_page_n_with_type_and_click(self, n):
            self.br.click("""//input[@class='pageNumBoxes']""")
            """http://apps.isiknowledge.com/summary.do?product=UA&search_mode=CitingArticles&qid=2&SID=3AhJNl1l543nfLJ2emf&page=2"""
            self.br.type("""//input[@class='pageNumBoxes']""", str(n) )
            self.br.click( """//input[@alt='Go to the page']""")

        def calc_num_pages(self, numCiting):
            numPages = numCiting / 50
            leftOver = numCiting % 50
            if leftOver > 0: numPages += 1
            return numPages
            

        def setup_save_options(self):
                # Setup save options		
		self.br.click( """//input[@value='allrecord']""" )  ## Choose to download "All items on page"
		self.br.click("abstract")                           ## *DE*SELECTS the "abstract" option
		self.br.select( "save_options" , "Save to BibTeX" ) ## BibTex format is easy to parse, so use it

######################################################################### end Class

def getPubMedPyData( baseDir, limit=None ):
        i=0
	data = defaultdict(dict)
	#items = ['AU','TI']
        items = ['TI']
	
	conv = { 'AU' : 'author', 'TI' : 'title' }
	files = os.listdir( os.sep.join( [ baseDir ,'TI'] ) )

	files.sort()

        if limit is None:   LIMIT = len(files)
        else:               LIMIT = limit

        print >>stderr, 'Reading PubMedPy Data...'
	
	for pmid in files[:LIMIT]:
                i += 1
                print >>stderr, i, 
		for item in items:
                        fileHandle = open( os.path.join( baseDir, item, pmid  ) ,'r' )
			fileData = fileHandle.read()
			data[pmid][ conv[item] ] = fileData
			fileHandle.close()
		data[pmid]['pmid'] = pmid


        print >>stderr, '\nFinished Reading PubMedPy Data...'
	
	return data


def get_done():
    done = [ x.partition('.')[0] for x in os.listdir( C.DATA_DIR ) ]    # Done pmids. Filenames have '.bib' extension, so must strip that off

    done.sort()
    
    try:
        last_one = done[-1]
    except:
        last_one = 'There is none to remove'
        
    done = done[:-1]    # Chop off the last one, it may not have been finished properly, so dont count it as done
    DONE = dict( [(x,True) for x in done] )  # Create a dictionary of completed PMIDs, for fast searching by PMID
    return DONE, last_one


#####################################################################################################################################
### DEBUGGING AND TESTING RELATED THINGS ###
############################################

def in_test(): return bool( 0 )

def debug():   return bool( 1 )

def print_debug(s):
    if debug():
        stderr.write(s)

#####################################################################################################################################

def main():
        browser = Browser()

        ### IF IN_TEST:    ###
        if in_test():
            browser.run_article_TEST('AMP-activated protein kinase phosphorylation of endothelial NO synthase.', "555")
            exit()
        ### END IF IN_TEST ###

        DONE, last = get_done()

        print >>stderr, '\nThese PMIDs have already been done:'
        print >>stderr, sorted( DONE.keys() )
        print >>stderr, '\nPrevious Last PMID (considered incomplete): \'' + last + "\'"
        print >>stderr, ''
        
	pubmed_data = getPubMedPyData( C.PUBMED_DATA_LOCATION, C.DATA_LIMIT )
	if C.DATA_LIMIT is None:    LAST = len(pubmed_data)
	else:			    LAST = C.DATA_LIMIT
	
	for pmid in sorted(pubmed_data.keys())[:LAST]:
                if pmid in DONE: continue
                article = pubmed_data[pmid]
		if debug():
                        print '\n'
                        print 'Next PubMed Article:', article
                browser.run_article( article )
        browser.stop()
        
#####################################################################################################################################

if __name__ == "__main__":
        pass
	main()
	
#####################################################################################################################################