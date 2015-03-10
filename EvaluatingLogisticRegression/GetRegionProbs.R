library(reshape2)
library(reshape)

# RScript GetRegionProbs.R '/Users/cmelton/Documents/Aptana Studio 3 Workspace/RegStatsCodeForPaper/EvaluatingLogisticRegression/TestData/Chrom22MutCovariates.tsv' '/Users/cmelton/Documents/Aptana Studio 3 Workspace/RegStatsCodeForPaper/EvaluatingLogisticRegression/TestData/newLRFit' '/Users/cmelton/Documents/Aptana Studio 3 Workspace/RegStatsCodeForPaper/EvaluatingLogisticRegression/TestData/Chrom22MutCovariates_Prob.tsv'

args<-commandArgs(TRUE)
inputFile<-args[1]
print(inputFile)
lrModelFilename<-args[2]
print(lrModelFilename)
outputFile<-args[3]
print(outputFile)
load(file=lrModelFilename)

getRegionProb<-function(filename, pid, lrFit, pids) {
  WindowCovariate<-read.delim(filename,header=T, sep="\t", stringsAsFactors=F)
  WindowCovariate$regionID<-paste(WindowCovariate$chrom, WindowCovariate$start, sep="_")
  
  # convert bpType to refPair
  WindowCovariate$refPair[WindowCovariate$bpType==1]<-"AT"
  WindowCovariate$refPair[WindowCovariate$bpType==2]<-"CG"
  WindowCovariate$refPair<-factor(WindowCovariate$refPair, levels<-c("AT", "CG"))
  
  # convert transcriptType
  WindowCovariate$transcriptType[WindowCovariate$transcriptType==0]<-"none"
  WindowCovariate$transcriptType[WindowCovariate$transcriptType==1]<-"coding_exon"
  WindowCovariate$transcriptType[WindowCovariate$transcriptType==2]<-"coding_intron"
  WindowCovariate$transcriptType[WindowCovariate$transcriptType==3]<-"noncoding_exon"
  WindowCovariate$transcriptType[WindowCovariate$transcriptType==4]<-"noncoding_intron"
  WindowCovariate$transcriptType<-factor(WindowCovariate$transcriptType, levels<-c("coding_exon", "coding_intron", "noncoding_exon", "noncoding_intron", "none"))
  
  # prob from glm
  WindowCovariate$logProb<-0
  WindowCovariate$pid<-factor(pid)
  WindowCovariate$pidProb<-predict.glm(object=lrFit, newdata=WindowCovariate, type="response")
  WindowCovariate$logProb<-log10(1-WindowCovariate$pidProb)
  
  # sum region log probs
  WindowCovariateMelted<-melt(WindowCovariate, id.vars=c("regionID"), measure.vars=c("logProb"))
  WindowCovariateCasted<-cast(WindowCovariateMelted, regionID~variable, fun.aggregate=sum)
  WindowCovariateCasted$prob<-(1-10^(WindowCovariateCasted$logProb))
  return(WindowCovariateCasted[, c("regionID", "prob")])
}

pids<-unique(as.character(newLRFit$data$pid))

regionProb<-getRegionProb(inputFile,pids[1],newLRFit, pids)
colnames(regionProb)<-c("regionID", pids[1])
for (pid in pids[2:length(pids)]) {
  newRegionProb<-getRegionProb(inputFile,pid,newLRFit, pids)
  colnames(newRegionProb)<-c("regionID", pid)
  regionProb<-merge(regionProb, newRegionProb, by="regionID")
}

write.csv(regionProb, file=outputFile, quote=F, row.names=F)
