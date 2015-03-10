import thread, datetime

# This function allows multiple threads to write to the same logfile, it prevents
# multiple threads from writing at the same time by adding a lock
class LogFile():
    def __init__(self, fileName):
        self.fileName=fileName
        self.lock=thread.allocate_lock()
        try:
            f = open(self.fileName, 'w')
            f.write("logfile for Cancer Pipeline")
            f.close()
        except:
            print "something is wrong with log file"
        
    def write(self, textToWrite):
        self.lock.acquire()
        try:
            f = open(self.fileName, 'a')
            # write time then add text to write
            towrite=datetime.datetime.now().strftime("%Y-%m-%d %H:%M")+"\t"+textToWrite
            print towrite
            f.write("\n"+towrite)
            if f: f.close()
        except:
            print "something is wrong with log file writing"
        self.lock.release()