from GridEngineJobThread_v2 import Job

# this is a class for running a single job in the pipeline, we would write many of these,
# one for each distinct job/component of the pipeline, it is subclassed off of Job which
# has some nice functions that allow you to write a shell script interpretable by the grid engine,
# each of us could write one of these for each component of the pipeline we are responsible for
class pipelineJob(Job):
    def __init__(self, grid, jobInfoDict):
        self.grid=grid
        self.jobInfoDict=jobInfoDict
        self.dependencies=[]
        super(pipelineJob, self).__init__(self.grid, self.jobInfoDict["scriptName"])
        self.started = False
        self.ready = False
        self.finished = False
        self.failed = False
        self.status = None
    
    # update job status
    def updateStatus(self, status):
        self.status = status
        self.grid.logFileWriter.write("status of "+self.name+" is "+ self.status)
        if status == "done":
            self.finished = True
            self.grid.logFileWriter.write("updated status of "+self.name+" to "+self.status)
        if status=="failed":
            self.failed = True
            self.grid.logFileWriter.write("updated status of "+self.name+" to "+self.status)
        
    
    # check to see that dependencies haven't failed
    def dependenciesOkay(self, JobsDict):
        for d in self.dependencies:
            if JobsDict[d].failed: return False
        return True
    
    # check if dependencies are finished
    def readyToRun(self, JobsDict):
        if self.ready: return True
        if not self.ready:
            for d in self.dependencies:
                if d in JobsDict:
                    if not JobsDict[d].finished or JobsDict[d].failed:
                        return False
        self.ready = True 
        return True
            
    # set dependency names
    def setDependencies(self, dependencyDict):
        if self.jobInfoDict["dependencies"]=="": return
        for dep in self.jobInfoDict["dependencies"].split("|"):
            depName=("pat_"+self.jobInfoDict["patientID"]+"_"+dep)
            if depName in dependencyDict:
                self.dependencies=self.dependencies+dependencyDict[depName]
        print self.name+": "+", ".join(self.dependencies)
    
    # start the job
    def start(self):
        self.started = True
        # write a shell script to run
        self.writeShellScript(self.jobInfoDict["scriptPath"], self.jobInfoDict["scriptName"], 
                              self.jobInfoDict["scriptTime"], self.jobInfoDict["scriptErrorFileDirectory"],
                              self.jobInfoDict["scriptOutputFileDirectory"], self.jobInfoDict["scriptCustomizations"],
                              self.jobInfoDict["scriptCommand"], self.jobInfoDict["scriptMemory"], 
                              self.jobInfoDict["scriptEmailAddress"])
        # create job
        self.createJob(self.jobInfoDict["scriptPath"], self.jobInfoDict["inputs"].split("|"))
        
        # run job
        self.runJob()