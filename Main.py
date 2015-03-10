'''
Created on Feb 25, 2015

@author: cmelton
'''

# This package runs recurrent mutation analysis from mutation and coverage file inputs.
# The inputs for mutations should be in BED format and the coverage files in WIG format.

# Imports
from optparse import OptionParser
import os, threading, subprocess, Queue
# from GridEngine.Job_v2 import pipelineJob as GridEngineJob 
# from GridEngine.GridEngine import GridEngine
import time

# GLOBAL VARIABLES
CONVERT_TO_WIG_PROGRAM="./Preprocessing/MergedMutToWig.py"
MUTATION_COVARIATE_PROGRAM="./CovariateStats/CovariatesForEveryPositionInWig.py"
WHOLEGENOME_COVARIATE_PROGRAM="./CovariateStats/WigStats_forSingleOverlapWig.py"
FORMAT_PATIENT_MODEL_DATA_PROGRAM="./CovariateStats/CombineMutAndGenomeCovariates.R"
COLLAPSE_MUT_COVARIATE_DATA_PROGRAM="./CovariateStats/MutCovariatesSummary.R"
LOGISTICREGRESSION_MODEL_PROGRAM="./CovariateStats/GenerateModel.R"
POISSONBINOMIAL_RECURRENT_PROB_PROGRAM="./PoissonBinomial/PoiBinProbsNew.R"
MERGE_MUTATIONS_PROGRAM="./MergeMutations.py"
COMPUTE_MERGED_MUTATION_PROBS_PROGRAM="./PoissonBinomial/GetRegionProbs.R"
DEFAULTGRIDSCRIPTPATH="/home/cmelton/"
DEFAULTGRIDOUTPUTPATH="/home/cmelton/"
DEFAULTGRIDERRORPATH="/home/cmelton/"
DEFAULTGRIDEMAILADDRESS="cmelton@stanford.edu"

# this function reads the input file and generates a dictionary of column 
# one as key and column two as value
def GetFiles(MutationFileListFile):
    f=open(MutationFileListFile, 'r')
    lines=f.read().split("\n")
    f.close()
    header=lines[0].strip().split("\t")
    lines=lines[1:]
    result={}
    for line in map(lambda x: x.strip().split("\t"), lines):
        result[line[0]]=dict(map(lambda i: (header[i], line[i]), range(len(line))))
    return result

class GridTools():
    def __init__(self):
        pass

    @staticmethod
    # determine if any jobs remain
    def JobsRemain(JobsDict):
        for jobName in JobsDict:
            job=JobsDict[jobName]
            if not job.finished: return True
        return False
    
    @staticmethod
    # update jobs dict with all finished jobs
    def UpdateJobsDict(JobsDict):
        for jobName in JobsDict:
            job=JobsDict[jobName]
            job.updateStatus()
    
    @staticmethod
    # run all jobs until none remain
    def RunJobs(JobsDict, grid):
        # while jobs remain check if jobs are ready and start them
        while (GridTools.JobsRemain(JobsDict)):
            # update Jobs Dict with newly completed Jobs
            grid.updateJobDict(JobsDict)
            # start jobs that are ready
            for jobName in JobsDict:
                job=JobsDict[jobName]
                if job.readyToRun(JobsDict) and not job.started:
                    job.start()
                    print job.name, job.started
            # wait before checking on jobs to run
            time.sleep(60)
    
    @staticmethod
    def CreateJobInfoDict(name, command):
        return {"dependencies": "",
                "scriptPath": DEFAULTGRIDSCRIPTPATH,
                "scriptName": name,
                "scriptTime": "24:00:00",
                "scriptErrorFileDirectory": DEFAULTGRIDOUTPUTPATH,
                "scriptOutputFileDirectory": DEFAULTGRIDERRORPATH,
                "scriptCustomizations": "",
                "scriptCommand": command,
                "scriptMemory": 6, 
                "scriptEmailAddress": DEFAULTGRIDEMAILADDRESS,
                "inputs": ""}

# a class to run a command through subprocess in a thread
class ThreadCommand(threading.Thread):
    """Simple class to run multiple commands in parallel"""
    def __init__(self, queue):
        threading.Thread.__init__(self)
        self.queue = queue

    def run(self):
        while True:
            command = self.queue.get()
            try: result=subprocess.check_output(command,shell=True, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                result = "\n".join(map(str, [e.cmd, e.returncode, e.output]))
                print "failed", result
            self.queue.task_done()

# runs multiple commands in parallel via the thread and subprocess modules
def RunCommands(commands, n, grid):
    if grid!=None:
        pass
#         JobsDict = dict(map(lambda x: (str(i), GridEngineJob(grid, GridTools.CreateJobInfoDict(str(i), commands[i]))), range(len(commands))))
#         GridTools.RunJobs(JobsDict, grid)
    else:
        queue = Queue.Queue()
        #spawn a pool of threads, and pass them queue instance 
        for i in range(min(len(commands), n)):
            t = ThreadCommand(queue)
            t.setDaemon(True)
            t.start()
        
        #populate queue with data   
        for command in commands:
            queue.put(command)
     
        #wait on the queue until everything has been processed     
        queue.join()    

# generates a combined covariate file to speed up processing
def CombineCovariates(CovariateFiles, CombinedCovariatePath):
    if not os.path.exists(CombinedCovariatePath):
        pass 

# for each mutation file get covariates for each mutation
def GetMutationCovariates(MutationFiles, CombinedCovariateFile, parallel, grid):
    # convert to wig files
    commands=[]
    for sample in MutationFiles.values():
        commands.append("python "+CONVERT_TO_WIG_PROGRAM+" --I "+sample["MutationFile"]+" --O "+sample["MutationWigFile"])
    RunCommands(commands, parallel, grid)
     
    # get covariates
    commands=[]
    for sample in MutationFiles.values():
        commands.append("python "+MUTATION_COVARIATE_PROGRAM+" --I "+sample["MutationWigFile"]+" --O "+sample["MutationCovariateFile"]+" --V "+CombinedCovariateFile)
    RunCommands(commands, parallel, grid)
    
    # combine into summary file
    commands=[]
    for sample in MutationFiles.values():
        commands.append("Rscript "+COLLAPSE_MUT_COVARIATE_DATA_PROGRAM+" "+sample["MutationCovariateFile"]+" "+sample["MutationCovariateSummaryFile"])
    RunCommands(commands, parallel, grid)
    
# get covariate totals for whole genome
def GetWholeGenomeCovariates(MutationFiles, CombinedCovariateFile, parallel, grid):
    commands=[]
    for sample in MutationFiles.values():
        commands.append("python "+WHOLEGENOME_COVARIATE_PROGRAM+" --I "+sample["CoverageWigFile"]+" --O "+sample["WGCovariateFile"]+" --V "+CombinedCovariateFile)
    RunCommands(commands, parallel, grid)
    
    # combine into summary file
    commands=[]
    for sample in MutationFiles.values():
        commands.append("Rscript "+FORMAT_PATIENT_MODEL_DATA_PROGRAM+" "+sample["WGCovariateFile"]+" "+sample["MutationCovariateSummaryFile"]+" "+sample["ModelData"]+" "+sample["pid"])
    RunCommands(commands, parallel, grid)
    
# fit the logistic regression sample specific probability model
def GenerateLRModel(MutationFiles, LRModelName):
    modeldata="|".join(map(lambda sample:sample["ModelData"], MutationFiles.values()))
    print modeldata
    subprocess.call("Rscript "+LOGISTICREGRESSION_MODEL_PROGRAM+" '"+modeldata+"' '"+LRModelName+"'", shell=True)

# combine mutations into a single file with chrom pos and count info
def MergeMutations(MutationFiles, MergedMutationFilename, splitByChrom=True):
    mutationFiles="||".join(map(lambda sample: sample["pid"]+"|"+sample["MutationFile"], MutationFiles.values()))
    if splitByChrom:
        subprocess.call("python "+MERGE_MUTATIONS_PROGRAM+" --I '"+mutationFiles+"' --O "+MergedMutationFilename+" --S T", shell=True)
    else:
        subprocess.call("python "+MERGE_MUTATIONS_PROGRAM+" --I '"+mutationFiles+"' --O "+MergedMutationFilename+" --S F", shell=True)

# Get Sample Specific Probabilities for each mutation
def GetSampleSpecificMutationProbs(MergedMutationFilename, LRModelName, parallel, grid):
    MergedMutationFilenames=["/".join(MergedMutationFilename.split("/")[0:-1])+"/"+x for x in os.listdir("/".join(MergedMutationFilename.split("/")[0:-1])) if MergedMutationFilename.split("/")[-1] in x]
    commands=[]
    for filename in MergedMutationFilenames:
        commands.append("Rscript "+COMPUTE_MERGED_MUTATION_PROBS_PROGRAM+" "+filename+" "+LRModelName+" "+filename+".prob.tsv")
    RunCommands(commands, parallel, grid)
        
# Combine sample specific probabilitites with Poisson Binomial to get Recurrence Probabilities
def ComputePoissonBinomialProbs(MergedMutationFilename, parallel, grid, regionSize="'1'"):
    MergedMutationProbFilenames=["/".join(MergedMutationFilename.split("/")[0:-1])+"/"+x for x in os.listdir("/".join(MergedMutationFilename.split("/")[0:-1])) if MergedMutationFilename.split("/")[-1] in x and ".prob.tsv" in x]
    commands=[]
    for filename in MergedMutationProbFilenames:
        numberOfMutationsPerSiteData=filename[0:-9]
        outputCSVFile="'"+filename+".poibin.csv"+"'"
        commands.append("Rscript "+POISSONBINOMIAL_RECURRENT_PROB_PROGRAM+" "+filename+" "+numberOfMutationsPerSiteData+" "+outputCSVFile+" "+regionSize)
    RunCommands(commands, parallel, grid)

# this functions gets the command line options for running the program
def getOptions():
    parser = OptionParser()
    parser.add_option("--M", dest = "MutationFileListFile", help = "this should be a tab delimited file with patient id and mutation file location",
                      metavar = "STRING", type = "string", default = "MutationFiles.tsv")
    parser.add_option("--C", dest = "CovariateFileListFile", help = "this should be a tab delimited file with covariate file name and file location",
                      metavar = "STRING", type = "string", default = "")
    parser.add_option("--CC", dest = "CombinedCovariateFile", help = "this should be a filename for the combined covariates",
                      metavar = "STRING", type = "string", default = "./CovariateStats/TestData/mergedBPRepTimingTranscript.wig.gz")#"./CovariateStats/TestData/mergedBPRepTimingTranscript.22.wig.gz")
    parser.add_option("--LR", dest = "LRModelName", help = "this should be a filename for the fitted logistic regression model",
                      metavar = "STRING", type = "string", default = "./Temp/LRModel")
    parser.add_option("--P", dest = "parallel", help = "the number of jobs to run in parallel",
                      metavar = "STRING", type = "string", default = "2")
    parser.add_option("--MF", dest = "MergedMutationFilename", help = "the path to the merged mutation file",
                      metavar = "STRING", type = "string", default = "./Temp/MergedMutations.tsv")
    parser.add_option("--G", dest = "grid", help = "T or F to specify use of the grid engine",
                      metavar = "STRING", type = "string", default = "F")
    parser.add_option("--L", dest = "logFilePath", help = "the path to a log file for the grid engine",
                      metavar = "STRING", type = "string", default = "")
    parser.add_option("--RS", dest = "regionSize", help = "",
                      metavar = "STRING", type = "string", default = "1")
    (options, args) = parser.parse_args()
    return options


if __name__ == '__main__':
    # get options
    options = getOptions()
    
    #initialize grid engine if needed
    grid=None
    if options.grid=="T":
        pass
#         grid = GridEngine(options.logFileFilePath)

    try:
        # get a dictionary of patient id and mutation file
        MutationFiles=GetFiles(options.MutationFileListFile)
        print MutationFiles
#         # Preprocess to get combined covariate file
#         CovariateFiles=GetFiles(options.CovariateFileListFile)
#         CombineCovariates(CovariateFiles, options.CombinedCovariateFile)
#          
#         # Get Mutation Covariates
#         GetMutationCovariates(MutationFiles, options.CombinedCovariateFile, int(options.parallel), grid)
#          
#         # Get Whole Genome Covariates
#         GetWholeGenomeCovariates(MutationFiles, options.CombinedCovariateFile, int(options.parallel), grid)
#          
#         # Generate LR Sample Specific Probability Model
#         GenerateLRModel(MutationFiles, options.LRModelName)
#          
#         # Get Merged Mutations
#         MergeMutations(MutationFiles, options.MergedMutationFilename, splitByChrom=True)
#          
#         # Get Sample Specific Probabilities
#         GetSampleSpecificMutationProbs(options.MergedMutationFilename, options.LRModelName, int(options.parallel), grid)
#          
        # Compute Poisson Binomial Recurrence Probabilities
        ComputePoissonBinomialProbs(options.MergedMutationFilename, int(options.parallel), grid, regionSize="'"+options.regionSize+"'")

    finally:
        # exit drmaa session
        if grid!=None:
            grid.exit()
    
        
    