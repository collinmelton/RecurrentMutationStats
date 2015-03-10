library(reshape)
library(reshape2)

args<-commandArgs(TRUE)
wg<-args[1] #"/Users/cmelton/Documents/AptanaStudio3Workspace/RecurrentMutationWrapper/Temp/Chrom22Covariates.tsv"#
mut<-args[2] #"/Users/cmelton/Documents/AptanaStudio3Workspace/RecurrentMutationWrapper/Temp/Chrom22MutCovariates.tsv.merged"#
output<-args[3] #"/Users/cmelton/Documents/AptanaStudio3Workspace/RecurrentMutationWrapper/Temp/ModelData"#
pid<-args[4] #"1" #

allPos<-read.delim(wg, header=T)
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
muts<-read.delim(mut, header=T)
muts$transcript_type<-"none"
muts$transcript_type[muts$coding_exon==1]<-"coding_exon"
muts$transcript_type[muts$coding_intron==1]<-"coding_intron"
muts$transcript_type[muts$noncoding_exon==1]<-"noncoding_exon"
muts$transcript_type[muts$noncoding_intron==1]<-"noncoding_intron"
muts$covered<-1
muts$refPair<-"NN"
muts$refPair[muts$bpType==1]<-"AT"
muts$refPair[muts$bpType==2]<-"CG"
formattedMuts<-muts[,c("refPair", "replicationTiming","transcript_type", "count")]

merged<-merge(formattedAllPos[, c("refPair", "replicationTiming","transcript_type", "count")], formattedMuts[, c("refPair", "replicationTiming","transcript_type", "count")], by=c("refPair", "replicationTiming","transcript_type"), all.x=T)
colnames(merged)[colnames(merged)=="count.x"]<-"all"
colnames(merged)[colnames(merged)=="count.y"]<-"mutant"
colnames(merged)[colnames(merged)=="transcript_type"]<-"transcriptType"
merged$mutant[is.na(merged$mutant)]<-0

merged$normal<-merged$all-merged$mutant
merged$refPair<-factor(merged$refPair)
merged$transcriptType<-factor(merged$transcriptType)

merged$prob<-merged$mutant/merged$all
merged$real<-T
merged$pid<-pid

merged$refPair<-as.character(merged$refPair)
merged$transcriptType<-as.character(merged$transcriptType)

write.table(merged, file = output, sep = "\t", quote = F, row.names = F, col.names = T)
