from optparse import OptionParser
import os

def run(inputFiles, outputFile):
    oldheader=""
    results={}
    for filename in inputFiles:
        f=open(filename, "r")
        header=f.readline()
        if filename != inputFiles[0]:
            if header!=oldheader: print "header error!!! headers don't match"
        data=f.read().split("\n")
        f.close()
        for line in data:
            vals=line.split("\t")
            if len(vals)>2:
                key=vals[0]
                count=int(vals[-1])
                middle="".join(map(lambda x: "\t"+x, vals[1:-1]))
                if key not in results:
                    results[key]={"middle": middle,
                                  "count": 0}
                results[key]["count"]+=count
        oldheader=header
    f=open(outputFile, "w")
    f.write(header)
    for key in results:
        f.write("\n"+key+results[key]["middle"]+"\t"+str(results[key]["count"]))
    f.close()
            
def getOptions():
    parser = OptionParser()
    parser.add_option("--P", dest="prefix", help = "path plus input file prefix",
                      metavar = "FILE", type = "string", default = "./TestData/CovariatesByChrom/")
    parser.add_option("--O", dest="outputFile", help="path and name of output file",
                      metavar="FILE", type="string", default="./TestData/MergedCovariates.txt")
    (options, args) = parser.parse_args()
    return options

# this is the main function
def runMain():
    options = getOptions()
    chroms=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]
    inputFiles=map(lambda x: os.path.join(options.prefix, "coverage.wig.corrected."+x+".NewWigStats.txt"), chroms)
    #inputFiles=[options.prefix, options.prefix]
    run(inputFiles, options.outputFile)    

runMain()