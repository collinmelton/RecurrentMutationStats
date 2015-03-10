import os
from optparse import OptionParser
import OverlapWig

class gencodeTranscriptWigs:
    def __init__(self, cExonFiles, cIntronFiles, ncExonFiles, ncIntronFiles):
        self.cExonFiles=cExonFiles
        self.cIntronFiles=cIntronFiles
        self.ncExonFiles=ncExonFiles
        self.ncIntronFiles=ncIntronFiles
        self.cExons=[]
        for f in self.cExonFiles:
            if os.path.exists(f):
                self.cExons.append(OverlapWig.OverlapWig(f))
        self.cIntrons=[]
        for f in self.cIntronFiles:
            if os.path.exists(f):
                self.cIntrons.append(OverlapWig.OverlapWig(f))
        self.ncExons=[]
        for f in self.ncExonFiles:
            if os.path.exists(f):
                self.ncExons.append(OverlapWig.OverlapWig(f))
        self.ncIntrons=[]
        for f in self.ncIntronFiles:
            if os.path.exists(f):
                self.ncIntrons.append(OverlapWig.OverlapWig(f))
    
    def getGencodeAnnotation(self, chrom, pos):
        for cExonWig in self.cExons:
            if cExonWig.getWigResults(chrom, pos)!=None:
                return(1)
        for cIntronWig in self.cIntrons:
            if cIntronWig.getWigResults(chrom, pos)!=None:
                return(2)
        for ncExonWig in self.ncExons:
            if ncExonWig.getWigResults(chrom, pos)!=None:
                return(3)
        for ncIntronWig in self.ncIntrons:
            if ncIntronWig.getWigResults(chrom, pos)!=None:
                return(4)
        return(0)

class fastaWig:
    def __init__(self, fastaWigFile):
        self.f=open(fastaWigFile, 'r')
        self.header=self.f.readline().strip()
        self.over=False
        print "header:", self.header
        self.chrom=self.header.split("chrom=")[1].split()[0]
        self.pos=int(self.header.split("start=")[1].split()[0])
        self.step=int(self.header.split("step=")[1].split()[0])
        
    def walk(self):
        if self.over: return None
        next=self.f.readline()
        if next=="": 
            self.over=True
            self.f.close()
            return None
        p= (next.strip(), self.chrom, self.pos)
        self.pos=self.pos+self.step
        #print p
        return p
    
    def getWigHeader(self):
        return self.header

class RepTimingWig:
    def __init__(self, filename, min, max, bins):
        #self.covWig=wig.Wig(filename)
        self.covWig=OverlapWig.OverlapWig(filename)
        self.min=min
        self.max=max
        self.bins=bins
        self.unit=(max-min)/bins
        
    def getBin(self, value):
        return int((float(value)-self.min)/self.unit)
    
    def getWigResults(self, chrom, pos):
        try:
            v = self.covWig.getWigResults(chrom, pos)
            if v==None: return -1
            b = self.getBin(v)
            return b
        except:
            print "Exception:"
            print (chrom, pos, self.covWig.getWigResults(chrom, pos))

def run(fastaWigFile, cExonFiles, cIntronFiles, ncExonFiles, ncIntronFiles, repTimingFile, repTimingMin, repTimingMax, repTimingNumBins, resultFile):
    print fastaWigFile
    print repTimingFile
    fw=fastaWig(fastaWigFile)
    print "coverage file in memory"
    rep=RepTimingWig(repTimingFile, repTimingMin, repTimingMax, repTimingNumBins)
    print "replication timing in memory"
    gencodeTranscripts=gencodeTranscriptWigs(cExonFiles, cIntronFiles, ncExonFiles, ncIntronFiles)
    
    f=open(resultFile, "w")
    print fw.getWigHeader()
    f.write(fw.getWigHeader())
    p=fw.walk()
    i=0
    while p!=None:
        i+=1
        val, chrom, pos=p
        if i%100000==0: print chrom, pos
        repTiming=str(rep.getWigResults(chrom, pos)) 
        #print chrom, pos, repTiming
        # get transcript annotations
        gencodeAnnotation=gencodeTranscripts.getGencodeAnnotation(chrom, pos)
        #print "_".join([str(val), str(repTiming), str(gencodeAnnotation)])
        f.write("\n"+"_".join([str(val), str(repTiming), str(gencodeAnnotation)]))
        p=fw.walk()
        #if i==100000: break
    f.close()
    

# this functions gets the command line options for running the program
def getOptions():
    parser = OptionParser()
    parser.add_option("--I", dest = "inputFile", help = "input file in wig format",
                      metavar = "FILE", type = "string", default = "/Users/cmelton/Downloads/hs37d5.CG_TA.1.wig")
    parser.add_option("--O", dest="outputFile", help = "output file is a tab delimited file of key and count",
                      metavar = "FILE", type = "string", default = "testout_overlapwig.out")
    parser.add_option("--D", dest = "WigDirectory", help = "path to directory of reference wig files"+
                      "concurrently", metavar = "FILE", default = "/Users/cmelton/Documents/Lab/SnyderLab/LocalMutationRate/WigFiles", type = "string")
    parser.add_option("--C", dest = "Chrom", help = "chromosome to use", metavar = "FILE", 
                  default = "1", type = "string")
    (options, args) = parser.parse_args()
    return options

# this is the main function
def runMain():
    # get options and defaults
    print "starting run"
    options = getOptions()
    directory=options.WigDirectory
    ncExonfilePrefixes=["3prime_overlapping_ncrna_exonInfo","IG_C_gene_exonInfo", "IG_D_gene_exonInfo",
            "IG_J_gene_exonInfo", "IG_V_gene_exonInfo","Mt_rRNA_exonInfo",
            "Mt_tRNA_exonInfo","TR_C_gene_exonInfo","TR_D_gene_exonInfo",
            "TR_J_gene_exonInfo","TR_V_gene_exonInfo","antisense_exonInfo",
            "lincRNA_exonInfo","miRNA_exonInfo","misc_RNA_exonInfo",
            "polymorphic_pseudogene_exonInfo","processed_transcript_exonInfo",
            "pseudogene_exonInfo",
            "rRNA_exonInfo","sense_intronic_exonInfo",
            "sense_overlapping_exonInfo","snRNA_exonInfo",
            "snoRNA_exonInfo", "noncodingtest"]
    ncIntronfilePrefixes=["3prime_overlapping_ncrna_intronInfo", 
            "IG_C_gene_intronInfo", "IG_V_gene_intronInfo", "TR_C_gene_intronInfo", 
            "TR_V_gene_intronInfo", "sense_intronic_intronInfo", "sense_overlapping_intronInfo", 
            "antisense_intronInfo", "lincRNA_intronInfo", "polymorphic_pseudogene_intronInfo", 
            "processed_transcript_intronInfo", "pseudogene_intronInfo", "noncodingintrontest"]
    cExonfilePrefixes=["protein_coding_exonInfo", "codingtest"]
    cIntronfilePrefixes=["protein_coding_intronInfo", "codingintrontest"]
    #gencodeFiles=map(lambda x: os.path.join(directory, x+"."+options.Chrom+".wig"), filePrefixes)
    cExonFiles=map(lambda x: os.path.join(directory, x+"."+options.Chrom+".wig"), cExonfilePrefixes)
    cIntronFiles=map(lambda x: os.path.join(directory, x+"."+options.Chrom+".wig"), cIntronfilePrefixes)
    ncExonFiles=map(lambda x: os.path.join(directory, x+"."+options.Chrom+".wig"), ncExonfilePrefixes)
    ncIntronFiles=map(lambda x: os.path.join(directory, x+"."+options.Chrom+".wig"), ncIntronfilePrefixes)
    
    print "about to run"
    run(options.inputFile, cExonFiles, cIntronFiles, ncExonFiles, ncIntronFiles, os.path.join(directory, "summary."+options.Chrom+".wig"), -4, 83, 87, options.outputFile)
    
runMain()
    
    
    
#print gencodeFiles




