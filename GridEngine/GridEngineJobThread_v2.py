from GridEngine import GridEngine

# This is the common super class to create and start jobs
class Job(object):
    def __init__(self, grid, name):
        self.name=name
        self.grid=grid
        self.drmaaSession=grid.drmaaSession
        #self.logFileWriter = self.grid.logFileWriter
        self.jobTemplate=None
        
    # write to the common log file
    def writeToLogFile(self, textToWrite):
        self.logFileWriter.write(textToWrite)
    
    # this function helps write a generic shell script for qsub
    def writeShellScript(self, filePath, scriptName, runTime, error, output, customizations, executionLine, memory, emailAddress):
        newline = "\n"
        f = open(filePath, "w")
        f.write("#!/bin/sh" + newline)
        # set the name of the job
        f.write("#$ -N " + scriptName + newline)
        # set max memory usage per slot
        f.write("#$ -l h_vmem=" + str(memory) + "G" + newline)
        # set max run time
        f.write("#$ -l h_rt="+runTime+ newline)    
        # send mail when job ends or aborts
        f.write("#$ -m ea" + newline)
        # specify email address
        f.write("#$ -M " + emailAddress + newline)
        # check for errors in job submission
        f.write("#$ -w e" + newline)
        # specify output and error files= names
        f.write("#$ -o " + output + newline)
        f.write("#$ -e " + error + newline)
        # add custom qsub options
        f.write(customizations + newline)
        # add code to be run
        f.write(executionLine)
        f.close()     
    
    # this method generates a jobtemplate object to be run
    # script path is the path to the script .sh .py or whatever, that is to be run, 
    # args is a list of arguments for the script
    def createJob(self, scriptPath, args):
        self.jobTemplate = self.grid.drmaaSession.createJobTemplate()
        self.jobTemplate.remoteCommand = scriptPath
        self.jobTemplate.args = args
        # the -b no option tells sge to not read the command as a binary file,
        # this means all the sge commented out options will be read, default is -b yes
        self.jobTemplate.nativeSpecification = "-b no"
        
    # this method rus the jobtemplate object
    def runJob(self):
        self.grid.startJob(self.jobTemplate, self.name)

    # this method marks the jobevent complete
    def cleanup(self):
        self.drmaaSession.deleteJobTemplate(self.jobTemplate)