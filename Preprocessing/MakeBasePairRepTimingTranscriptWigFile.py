import os, SourceWig
from optparse import OptionParser
import OverlapWig, gzip

class epigenomeWigs:
    def __init__(self, epigenomeFiles):
        self.epigenomeFiles=epigenomeFiles
        self.epigenomes=[]
        for f in self.epigenomeFiles:
            if os.path.exists(f):
                self.epigenomes.append(OverlapWig.OverlapWig(f))

    def getWigData(self, w, chrom, pos):
        result=w.getWigResults(chrom, pos)
        if result==None: result="-1"
        return str(result)

    def getAnnotation(self, chrom, pos):
        return "_".join(map(lambda x: self.getWigData(x, chrom, pos), self.epigenomes))

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
        if ".gz" in fastaWigFile:
            self.f=gzip.open(fastaWigFile, 'r')
        else:
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

def run(fastaWigFile, cExonFiles, cIntronFiles, ncExonFiles, ncIntronFiles, repTimingFile, 
        repTimingMin, repTimingMax, repTimingNumBins, resultFile, epigenomeFiles):
    print fastaWigFile
    print repTimingFile
    fw=fastaWig(fastaWigFile)
    print "fasta file loaded"
    rep=RepTimingWig(repTimingFile, repTimingMin, repTimingMax, repTimingNumBins)
    print "replication timing file loaded"
    gencodeTranscripts=gencodeTranscriptWigs(cExonFiles, cIntronFiles, ncExonFiles, ncIntronFiles)
    epigenomes=epigenomeWigs(epigenomeFiles)
    f=open(resultFile, "w")
    print fw.getWigHeader()
    f.write(fw.getWigHeader())
    p=fw.walk()
    i=0
#     print "starting test"
#     print gencodeTranscripts.getGencodeAnnotation("22", 25150831)
#     print gencodeTranscripts.getGencodeAnnotation("22", 25150832)
#     print gencodeTranscripts.getGencodeAnnotation("22", 25150833)
#     print gencodeTranscripts.getGencodeAnnotation("22", 25150831+81)
#     print gencodeTranscripts.getGencodeAnnotation("22", 25150831+80)
#     print gencodeTranscripts.getGencodeAnnotation("22", 25150831+82)
    
    while p!=None:
        i+=1
        val, chrom, pos=p
#         if i%100==0:print chrom, pos, val
#         if i%100000==0: 
#             break
        repTiming=str(rep.getWigResults(chrom, pos)) 
        #print chrom, pos, repTiming
        # get transcript annotations
        gencodeAnnotation=gencodeTranscripts.getGencodeAnnotation(chrom, pos)
        epiAnnotation=epigenomes.getAnnotation(chrom, pos)
        #print "_".join([str(val), str(repTiming), str(gencodeAnnotation)])
        if epiAnnotation!="":
            f.write("\n"+"_".join([str(val), str(repTiming), str(gencodeAnnotation), epiAnnotation]))
        else:
            f.write("\n"+"_".join([str(val), str(repTiming), str(gencodeAnnotation)]))
        p=fw.walk()
#        if i==100000: break
    f.close()
    

# this functions gets the command line options for running the program
def getOptions():
    parser = OptionParser()
    parser.add_option("--I", dest = "inputFile", help = "input file in wig format",
                      metavar = "FILE", type = "string", default = "./AnnotationWigs/humang1kv37.singlets22.wig.gz")#"/Users/cmelton/Downloads/hs37d5.CG_TA.1.wig")
    parser.add_option("--O", dest="outputFile", help = "output file is a tab delimited file of key and count",
                      metavar = "FILE", type = "string", default = "overlap.22.out")
    parser.add_option("--D", dest = "WigDirectory", help = "path to directory of reference wig files"+
                      "concurrently", metavar = "FILE", default = "./AnnotationWigs/", type = "string") #"/srv/gsfs0/projects/snyder/collinmelton/RecurrentMutationDecectionAlgorithm/Wigs/", type = "string" )#
    parser.add_option("--C", dest = "Chrom", help = "chromosome to use", metavar = "FILE", 
                  default = "22", type = "string")
    parser.add_option("--E", dest = "useEpi", help = "whether to use epigenome wigs", metavar = "FILE", 
                  default = "T", type = "string")
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
    cExonFiles=map(lambda x: os.path.join(directory+"", x+"."+options.Chrom+".wig.gz"), cExonfilePrefixes)
    cIntronFiles=map(lambda x: os.path.join(directory+"", x+"."+options.Chrom+".wig.gz"), cIntronfilePrefixes)
    ncExonFiles=map(lambda x: os.path.join(directory+"", x+"."+options.Chrom+".wig.gz"), ncExonfilePrefixes)
    ncIntronFiles=map(lambda x: os.path.join(directory+"", x+"."+options.Chrom+".wig.gz"), ncIntronfilePrefixes)
    
    # epigenome files
    if options.useEpi=="T":
        epigenomefilePrefixes=["E001_15_coreMarks_mnemonics.bed","E002_15_coreMarks_mnemonics.bed","E003_15_coreMarks_mnemonics.bed","E004_15_coreMarks_mnemonics.bed","E005_15_coreMarks_mnemonics.bed","E006_15_coreMarks_mnemonics.bed","E007_15_coreMarks_mnemonics.bed","E008_15_coreMarks_mnemonics.bed","E009_15_coreMarks_mnemonics.bed","E010_15_coreMarks_mnemonics.bed","E011_15_coreMarks_mnemonics.bed","E012_15_coreMarks_mnemonics.bed","E013_15_coreMarks_mnemonics.bed","E014_15_coreMarks_mnemonics.bed","E015_15_coreMarks_mnemonics.bed","E016_15_coreMarks_mnemonics.bed","E017_15_coreMarks_mnemonics.bed","E018_15_coreMarks_mnemonics.bed","E019_15_coreMarks_mnemonics.bed","E020_15_coreMarks_mnemonics.bed","E021_15_coreMarks_mnemonics.bed","E022_15_coreMarks_mnemonics.bed","E023_15_coreMarks_mnemonics.bed","E024_15_coreMarks_mnemonics.bed","E025_15_coreMarks_mnemonics.bed","E026_15_coreMarks_mnemonics.bed","E027_15_coreMarks_mnemonics.bed","E028_15_coreMarks_mnemonics.bed","E029_15_coreMarks_mnemonics.bed","E030_15_coreMarks_mnemonics.bed","E031_15_coreMarks_mnemonics.bed","E032_15_coreMarks_mnemonics.bed","E033_15_coreMarks_mnemonics.bed","E034_15_coreMarks_mnemonics.bed","E035_15_coreMarks_mnemonics.bed","E036_15_coreMarks_mnemonics.bed","E037_15_coreMarks_mnemonics.bed","E038_15_coreMarks_mnemonics.bed","E039_15_coreMarks_mnemonics.bed","E040_15_coreMarks_mnemonics.bed","E041_15_coreMarks_mnemonics.bed","E042_15_coreMarks_mnemonics.bed","E043_15_coreMarks_mnemonics.bed","E044_15_coreMarks_mnemonics.bed","E045_15_coreMarks_mnemonics.bed","E046_15_coreMarks_mnemonics.bed","E047_15_coreMarks_mnemonics.bed","E048_15_coreMarks_mnemonics.bed","E049_15_coreMarks_mnemonics.bed","E050_15_coreMarks_mnemonics.bed","E051_15_coreMarks_mnemonics.bed","E052_15_coreMarks_mnemonics.bed","E053_15_coreMarks_mnemonics.bed","E054_15_coreMarks_mnemonics.bed","E055_15_coreMarks_mnemonics.bed","E056_15_coreMarks_mnemonics.bed","E057_15_coreMarks_mnemonics.bed","E058_15_coreMarks_mnemonics.bed","E059_15_coreMarks_mnemonics.bed","E061_15_coreMarks_mnemonics.bed","E062_15_coreMarks_mnemonics.bed","E063_15_coreMarks_mnemonics.bed","E065_15_coreMarks_mnemonics.bed","E066_15_coreMarks_mnemonics.bed","E067_15_coreMarks_mnemonics.bed","E068_15_coreMarks_mnemonics.bed","E069_15_coreMarks_mnemonics.bed","E070_15_coreMarks_mnemonics.bed","E071_15_coreMarks_mnemonics.bed","E072_15_coreMarks_mnemonics.bed","E073_15_coreMarks_mnemonics.bed","E074_15_coreMarks_mnemonics.bed","E075_15_coreMarks_mnemonics.bed","E076_15_coreMarks_mnemonics.bed","E077_15_coreMarks_mnemonics.bed","E078_15_coreMarks_mnemonics.bed","E079_15_coreMarks_mnemonics.bed","E080_15_coreMarks_mnemonics.bed","E081_15_coreMarks_mnemonics.bed","E082_15_coreMarks_mnemonics.bed","E083_15_coreMarks_mnemonics.bed","E084_15_coreMarks_mnemonics.bed","E085_15_coreMarks_mnemonics.bed","E086_15_coreMarks_mnemonics.bed","E087_15_coreMarks_mnemonics.bed","E088_15_coreMarks_mnemonics.bed","E089_15_coreMarks_mnemonics.bed","E090_15_coreMarks_mnemonics.bed","E091_15_coreMarks_mnemonics.bed","E092_15_coreMarks_mnemonics.bed","E093_15_coreMarks_mnemonics.bed","E094_15_coreMarks_mnemonics.bed","E095_15_coreMarks_mnemonics.bed","E096_15_coreMarks_mnemonics.bed","E097_15_coreMarks_mnemonics.bed","E098_15_coreMarks_mnemonics.bed","E099_15_coreMarks_mnemonics.bed","E100_15_coreMarks_mnemonics.bed","E101_15_coreMarks_mnemonics.bed","E102_15_coreMarks_mnemonics.bed","E103_15_coreMarks_mnemonics.bed","E104_15_coreMarks_mnemonics.bed","E105_15_coreMarks_mnemonics.bed","E106_15_coreMarks_mnemonics.bed","E107_15_coreMarks_mnemonics.bed","E108_15_coreMarks_mnemonics.bed","E109_15_coreMarks_mnemonics.bed","E110_15_coreMarks_mnemonics.bed","E111_15_coreMarks_mnemonics.bed","E112_15_coreMarks_mnemonics.bed","E113_15_coreMarks_mnemonics.bed","E114_15_coreMarks_mnemonics.bed","E115_15_coreMarks_mnemonics.bed","E116_15_coreMarks_mnemonics.bed","E117_15_coreMarks_mnemonics.bed","E118_15_coreMarks_mnemonics.bed","E119_15_coreMarks_mnemonics.bed","E120_15_coreMarks_mnemonics.bed","E121_15_coreMarks_mnemonics.bed","E122_15_coreMarks_mnemonics.bed","E123_15_coreMarks_mnemonics.bed","E124_15_coreMarks_mnemonics.bed","E125_15_coreMarks_mnemonics.bed","E126_15_coreMarks_mnemonics.bed","E127_15_coreMarks_mnemonics.bed","E128_15_coreMarks_mnemonics.bed","E129_15_coreMarks_mnemonics.bed"]
        epigenomeFiles=map(lambda x: os.path.join(directory+"EpigenomeWigs/", x+"."+options.Chrom+".wig.gz"), epigenomefilePrefixes)
    else:
        epigenomeFiles=[]
    print "about to run"
    run(options.inputFile, cExonFiles, cIntronFiles, ncExonFiles, ncIntronFiles, os.path.join(directory, "summary."+options.Chrom+".wig.gz"), -4, 83, 87, options.outputFile, epigenomeFiles)
    
runMain()
    
