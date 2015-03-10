# import WigTools.SourceWig as SourceWig
# import WigTools.OverlapWig as OverlapWig
import SourceWig, OverlapWig
from optparse import OptionParser

class Coverage:
    def __init__(self, coverageFile):
        self.cov=SourceWig.SourceWig(coverageFile)
        self.p=None
        self.header=self.cov.currentHeader
#        self.oldHeader=None
        
    def walk(self):
        self.p=self.cov.walk()
        newHeader=False
        if self.cov.currentHeader!=self.header:
            newHeader=True
#        self.oldHeader=self.header
        self.header=self.cov.currentHeader
        return newHeader, self.p
    
    def header_to_string(self):
        return("fixedStep chrom="+ str(self.cov.currentHeader.chrom)+" start="+str(self.cov.currentHeader.start)+" step="+str(self.cov.currentHeader.step))

def run(coverageFile, overlapWigFile, resultFile):
    print coverageFile
    print overlapWigFile
    cov=Coverage(coverageFile)
    overlap=OverlapWig.OverlapWig(overlapWigFile)
    
    newHeader, p=cov.walk()
    i=0
    f=open(resultFile, "w")
    f.write("chrom\tstart\tcovered\tbpType\treplicationTiming\ttranscriptType")
    while p!=None:
#         if newHeader:
#             overlap.currentPosition=overlap.currentHeader.start
        i+=1
        val, chrom, pos=p
        if i%100000==0: print chrom, pos 
        annotation=overlap.getWigResults(chrom, pos)
        if annotation==None:
            print "None Annotation:", cov.header.start, chrom, pos, overlap.currentHeader.chrom,overlap.currentHeader.start, overlap.nextHeader.chrom, overlap.nextHeader.start 
            annotation="NA_NA_NA"
        #print chrom, pos, annotation
        key=str(val)+"_"+annotation
        if "\x00" not in key:
            #if newHeader:
            if f.tell()!=0:f.write("\n")
            #    f.write(cov.header_to_string())
            f.write(str(cov.header.chrom)+"\t"+str(cov.header.start)+"\t"+"\t".join(key.split("_")))            
        else:
            print "weirdo key:", chrom, pos, val, key, str(val), str(key)
        newHeader, p=cov.walk()
    f.close()
    

# this functions gets the command line options for running the program
def getOptions():
    parser = OptionParser()
    parser.add_option("--I", dest = "inputFile", help = "input file in wig format",
                      metavar = "FILE", type = "string", default = "./TestData/WholeGenomeWindows_cds_chr22_10.txt.wig")
    parser.add_option("--O", dest="outputFile", help = "output file is a tab delimited file of key and count",
                      metavar = "FILE", type = "string", default = "./TestData/Chrom22CDSMutCovariates.tsv")
    parser.add_option("--V", dest = "overlapWigFile", help = "path to overlap wig file"+
                      "concurrently", metavar = "FILE", default = "./TestData/mergedBPRepTimingTranscript.22.wig", type = "string")
    (options, args) = parser.parse_args()
    return options

# this is the main function
def runMain():
    # get options and defaults
    print "starting run"
    options = getOptions()
    print "about to run"
    print options.inputFile, options.overlapWigFile, options.outputFile
    run(options.inputFile, options.overlapWigFile, options.outputFile)
    
runMain()
    
    
    
#print gencodeFiles




