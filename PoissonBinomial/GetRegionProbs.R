library(reshape2)
library(reshape)
library(parallel)

args<-commandArgs(TRUE)
inputFile<-args[1] #"/Users/cmelton/Documents/AptanaStudio3Workspace/RecurrentMutationWrapper/Temp/MergedMutations.tsv_22.tsv" #
print(inputFile)
lrModelFilename<-args[2] #"/Users/cmelton/Documents/AptanaStudio3Workspace/RecurrentMutationWrapper/Temp/LRModel" #
print(lrModelFilename)
outputFile<-args[3]
print(outputFile)
load(file=lrModelFilename)

getRegionProb<-function(filename, pid, fittedmodel) {
  WindowCovariate<-read.delim(filename,header=T, sep="\t", stringsAsFactors=F)
  WindowCovariate$regionID<-paste(WindowCovariate$chrom, WindowCovariate$pos, sep="_")
  
  # convert bpType to refPair
  WindowCovariate$refPair<-"NN"
  WindowCovariate$refPair[WindowCovariate$bpType==1]<-"AT"
  WindowCovariate$refPair[WindowCovariate$bpType==2]<-"CG"
#  WindowCovariate$refPair<-factor(WindowCovariate$refPair, levels<-c("AT", "CG"))
  WindowCovariate<-WindowCovariate[WindowCovariate$refPair!="NN",]
  
  # convert transcriptType
  WindowCovariate$transcriptType[WindowCovariate$transcriptType==0]<-"none"
  WindowCovariate$transcriptType[WindowCovariate$transcriptType==1]<-"coding_exon"
  WindowCovariate$transcriptType[WindowCovariate$transcriptType==2]<-"coding_intron"
  WindowCovariate$transcriptType[WindowCovariate$transcriptType==3]<-"noncoding_exon"
  WindowCovariate$transcriptType[WindowCovariate$transcriptType==4]<-"noncoding_intron"
#  WindowCovariate$transcriptType<-factor(WindowCovariate$transcriptType, levels<-c("coding_exon", "coding_intron", "noncoding_exon", "noncoding_intron", "none"))
  
  # prob from glm
  WindowCovariate$logProb<-0
  WindowCovariate$pid<-pid #factor(pid)
  chunksize<-50000
  chunks<-lapply(1:ceiling(dim(WindowCovariate)[1]/chunksize), function(i) ((i-1)*chunksize+1):(min((i)*chunksize, dim(WindowCovariate)[1])))
  fittedprobs<-lapply(chunks, function(chunk) predict.glm(object=fittedmodel, newdata=WindowCovariate[chunk,], type="response"))
  WindowCovariate$pidProb<-Reduce(function(x,y) c(x,y), fittedprobs)
  WindowCovariate$logProb<-log10(1-WindowCovariate$pidProb)
  
  # sum region log probs
  WindowCovariateMelted<-melt(WindowCovariate, id.vars=c("regionID"), measure.vars=c("logProb"))
  WindowCovariateCasted<-cast(WindowCovariateMelted, regionID~variable, fun.aggregate=sum)
  WindowCovariateCasted$prob<-(1-10^(WindowCovariateCasted$logProb))
  return(WindowCovariateCasted[, c("regionID", "prob")])
}

pids<-unique(fittedmodel$data$pid)

mergeCols<-function(regionProb, newRegionProb) {
  print(dim(newRegionProb))
  return(merge(regionProb, newRegionProb, by="regionID"))
}
setColName<-function(newRegionProb, pid) {
  colnames(newRegionProb)<-c("regionID", pid)
  return(newRegionProb)
}

fittedProbs<-mclapply(pids, function(pid) getRegionProb(inputFile,pids[1],fittedmodel), mc.cores = 4)
for (i in 1:length(fittedProbs)) {
  fittedProbs[[i]]<-setColName(fittedProbs[[i]], pids[i]) 
}
regionProb<-Reduce(mergeCols, fittedProbs)

write.csv(regionProb, file=outputFile, quote=F, row.names=F)
