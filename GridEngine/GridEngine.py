
import drmaa
import LogFile

# class drmaaClass():
#     def __init__(self):
#         self.name="fake"
#         
#     def Session(self):
#         return "fake session"

# drmaa=drmaaClass()

# this class handles the drmaa session and allows for control of starting and checking
# on grid engine jobs            
class GridEngine():
    def __init__(self, logFileFilePath):
        self.logFileWriter = LogFile.LogFile(logFileFilePath)
        self.drmaaSession = drmaa.Session()
        self.drmaaSession.initialize()
        self.logFileWriter.write("initializing drmaa session")
        self.currentJobs=set([])
        self.jobIDNameDict={}
    
    def exit(self):
        self.drmaaSession.exit()
        self.logFileWriter.write("exiting drmaa session")
    
    # start a job and add to current jobs
    def startJob(self, jobTemplate, jobName):
        jobid = self.drmaaSession.runJob(jobTemplate)
        self.currentJobs.add(jobid)
        self.jobIDNameDict[jobid]=jobName
        self.logFileWriter.write("starting job: "+jobName)
        
    # check status of current jobs and update job dictionary
    def updateJobDict(self, jobDict):
        self.logFileWriter.write("updating job status")
        self.logFileWriter.write("current jobs: "+",".join(map(str, list(self.currentJobs))))
        toremove=set([])
        for id in self.currentJobs:
            status=self.drmaaSession.jobStatus(id)
            jobDict[self.jobIDNameDict[id]].updateStatus(status)
            # remove from current jobs if done or failed
            if status=="done" or status=="failed":
                toremove.add(id)
        for id in toremove:
            self.currentJobs.remove(id)
            jobDict[self.jobIDNameDict[id]].cleanup()