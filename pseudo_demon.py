import time 
import sel_isi 

# start after the 90 minute delay
#time.sleep(25*60)


fname = 'crash_log.txt'

goon = True 
counter =0 
crash_time_old = time.time()


while goon== True:
    try:
        sel_isi.main()

    except:
        crash_time_new = time.time()
        elapsed_time = (crash_time_new - crash_time_old)
        counter +=1

        if (elapsed_time < 90) & (counter < 2 ):
            print "Crash occured in the last 90 sec - restart in 10 s"
            time.sleep(10)
        if (counter >=2):
            print "Crashed too many times, wait it out for 92 minutes"
            counter = 0 
            time.sleep(92*60)




