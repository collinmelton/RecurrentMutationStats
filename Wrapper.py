
import subprocess

# Notes: This wrapper is to show how the different modules fit together. In practice we 
# integrated these scripts in a grid engine environment so as to run jobs in parallel.

# Step 1: Get covariate (replication timing, transcript type, base pair composition) data
# for all individual mutations, all regions of interest for statistical testing, and 
# bulk statistics for the rest of the genome.

# 1a, covariate data for all individual mutations
inputFile= "./CovariateStats/TestData/01eef340-598c-4205-a990-cec190ac2ca5.unannotated.varscan.mutect.merge.wig_22.wig"# input file of genomic positions in wig format
outputFile="./CovariateStats/TestData/Chrom22MutCovariates.tsv" # output file in tsv format
covariateFile="./CovariateStats/TestData/mergedBPRepTimingTranscript.22.wig" # file with covariate data in wig format
subprocess.call("python ./CovariateStats/CovariatesForEveryPositionInWig.py --I "+inputFile+" --O "+outputFile+" --V "+covariateFile, shell=True)
 
# 1b, covariate data for all regions of interest for statistical testing
inputFile= "./CovariateStats/TestData/WholeGenomeWindows_1a.1b.1d.1f.2a.2b.2c.3a.3b.4.5_chr22_1.txt.wig"# input file of genomic positions in wig format
outputFile="./CovariateStats/TestData/Chrom22MutCovariates.tsv" # output file in tsv format
covariateFile="./CovariateStats/TestData/mergedBPRepTimingTranscript.22.wig" # file with covariate data in wig format
subprocess.call("python ./CovariateStats/CovariatesForEveryPositionInWig.py --I "+inputFile+" --O "+outputFile+" --V "+covariateFile, shell=True)
 
# 1c, covariate data bulk statistics for whole genome, the coverage file is an output of mutect
# which describes which sites have enough coverage to be called a mutation
inputFile= "./CovariateStats/TestData/coverage.wig.corrected.22.wig"# input file of genomic positions in wig format
outputFile="./CovariateStats/TestData/Chrom22Covariates.tsv" # output file in tsv format
covariateFile="./CovariateStats/TestData/mergedBPRepTimingTranscript.22.wig" # file with covariate data in wig format
subprocess.call("python ./CovariateStats/WigStats_forSingleOverlapWig.py --I "+inputFile+" --O "+outputFile+" --V "+covariateFile, shell=True)

# Step 2: Compute Logistic Regression model to get probability of mutation conditioned on covariates
print "step 2"
inputPathToFolderWithCovariatesForAllPositions="./LogisticRegression/TestData/MutCovariates/'" # outputs from 1c, files in folder must be named as follows: patientID.mergedWigResults
inputPathToFolderWithCovariatesForAllMutations="'./LogisticRegression/TestData/AllSiteCovariates/'" # outputs from 1a, files in folder must be named as follows: patientID.covariate_stats.txt
modelNamePath='./LogisticRegression/TestData/newLRFit' # name of file that will contain saved R model
tsvFileWithPatientIdentifiers="'./LogisticRegression/TestData/pids_test.tsv'" # file with names of patient ids used
subprocess.call("Rscript ./LogisticRegression/FormatDataAndComputeLRModel.R "+inputPathToFolderWithCovariatesForAllPositions+" "+inputPathToFolderWithCovariatesForAllMutations+" "+modelNamePath+" "+tsvFileWithPatientIdentifiers, shell=True)

# Step 3: Use logistic regression model to compute a table of sample specific conditional probabilities
print "step 3"
SitesToBeTestedCovariates="'./EvaluatingLogisticRegression/TestData/Chrom22MutCovariates.tsv'"
model="'./EvaluatingLogisticRegression/TestData/newLRFit'"
outputTSV="'./EvaluatingLogisticRegression/TestData/Chrom22MutCovariates_Prob.tsv'"
subprocess.call("RScript ./EvaluatingLogisticRegression/GetRegionProbs.R "+SitesToBeTestedCovariates+" "+model+" "+outputTSV, shell=True)

# Step 4: Calculate probability of k or more mutations per sites from step 3 using the Poisson Binomial 
print "step 4"
probabilityTSVFile="'./PoissonBinomial/TestData/Chrom22MutCovariates_Prob.tsv'"
numberOfMutationsPerSiteData="'./PoissonBinomial/TestData/WholeGenomeWindows_1a.1b.1d.1f.2a.2b.2c.3a.3b.4.5_chr22_1.txt'"
outputCSVFile="'./PoissonBinomial/TestData/PoiBinProbs.csv'"
regionSize="'1'"
subprocess.call("RScript ./PoissonBinomial/PoiBinProbs.R "+probabilityTSVFile+" "+numberOfMutationsPerSiteData+" "+outputCSVFile+" "+regionSize, shell=True)
