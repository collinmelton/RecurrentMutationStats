

# Rscript FormatDataAndComputeLRModel.R '~/Documents/Lab/SnyderLab/LocalMutationRate/PatientMergedWigStats/cleanstats/' '~/Documents/Lab/SnyderLab/LocalMutationRate/MutationWigFile/cleanstats_withRefPair/' '/Users/cmelton/Documents/Aptana Studio 3 Workspace/RegStatsCodeForPaper/LogisticRegression/TestData/newLRFit' '/Users/cmelton/Documents/Aptana Studio 3 Workspace/RegStatsCodeForPaper/LogisticRegression/TestData/pids_test.tsv'

args<-commandArgs(TRUE)
inputPathAllPos<-args[1] #"~/Documents/Lab/SnyderLab/LocalMutationRate/PatientMergedWigStats/cleanstats/"
print(inputPathAllPos)
inputPathMuts<-args[2] #"~/Documents/Lab/SnyderLab/LocalMutationRate/MutationWigFile/cleanstats_withRefPair/"
print(inputPathMuts)
lrModelFilename<-args[3] #"~/Documents/Lab/SnyderLab/LocalMutationRate/LogisticRegression/newLRFit"
print(lrModelFilename)
pidFile<-args[4] # "/Users/cmelton/Documents/Aptana Studio 3 Workspace/RegStatsCodeForPaper/LogisticRegression/pids.tsv"
print(pidFile)

pidDF<-read.delim(pidFile, header = T, sep = "\t", stringsAsFactor=F)
pids=pidDF$pid

# this function gets train data for a single pid
getMutNormalData<-function(pid, inputPathAllPos, inputPathMuts) {
  allPos<-read.delim(paste(inputPathAllPos, pid, ".mergedWigResults", sep=""), header=T)
  allPos$transcript_type<-"none"
  allPos$transcript_type[allPos$transcriptType==1]<-"coding_exon"
  allPos$transcript_type[allPos$transcriptType==2]<-"coding_intron"
  allPos$transcript_type[allPos$transcriptType==3]<-"noncoding_exon"
  allPos$transcript_type[allPos$transcriptType==4]<-"noncoding_intron"
  allPos$refPair<-"NN"
  allPos$refPair[allPos$bpType==1]<-"AT"
  allPos$refPair[allPos$bpType==2]<-"CG"
  formattedAllPos<-allPos[allPos$covered==1&allPos$bpType!=0,c("refPair", "replicationTiming","transcript_type", "count")]
  
  # process muts
  muts<-read.delim(paste(inputPathMuts, pid, ".covariate_stats.txt", sep=""), header=T)
  muts$transcript_type<-"none"
  muts$transcript_type[muts$coding_exon==1]<-"coding_exon"
  muts$transcript_type[muts$coding_intron==1]<-"coding_intron"
  muts$transcript_type[muts$noncoding_exon==1]<-"noncoding_exon"
  muts$transcript_type[muts$noncoding_intron==1]<-"noncoding_intron"
  muts$covered<-1
  formattedMuts<-muts[,c("refPair", "replicationTiming","transcript_type", "count")]
  
  merged<-merge(formattedAllPos[, c("refPair", "replicationTiming","transcript_type", "count")], formattedMuts[, c("refPair", "replicationTiming","transcript_type", "count")], by=c("refPair", "replicationTiming","transcript_type"), all.x=T)
  colnames(merged)[colnames(merged)=="count.x"]<-"all"
  colnames(merged)[colnames(merged)=="count.y"]<-"mutant"
  colnames(merged)[colnames(merged)=="transcript_type"]<-"transcriptType"
  merged$mutant[is.na(merged$mutant)]<-0
  # normal count is all - mutant
  merged$normal<-merged$all-merged$mutant
  merged$refPair<-factor(merged$refPair)
  merged$transcriptType<-factor(merged$transcriptType)
  
  merged$prob<-merged$mutant/merged$all
  merged$real<-T
  merged$pid<-pid
  return(merged)
}

# get merged data
allMerged<-do.call(rbind, lapply(pids, function(pid) getMutNormalData(pid, inputPathAllPos, inputPathMuts)))
allMerged$refPair<-as.character(allMerged$refPair)
allMerged$transcriptType<-as.character(allMerged$transcriptType)

# fit model
newLRFit<-glm(cbind(mutant, normal) ~ pid+refPair*pid+refPair+replicationTiming+pid*replicationTiming+replicationTiming+transcriptType, family = binomial, data=allMerged, x=F, y=F)

# save model
save(newLRFit, file=lrModelFilename)

