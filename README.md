# RecurrentMutationStats
This contains code for statistical analysis of recurrent mutations in whole genome sequencing data. 

# Contents
- [Overview](#overview)
- [Usage](#usage)
- [Steps](#steps) 

# Overview

This project includes code to perform the recurrent mutation analysis described in [Melton et al. Nature Genetics 2015](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4485503/)
This analysis is comprised of two steps: (1) Build a sample and genomic location specific mutation probability model and (2) use the Poisson Binomial to compute the probability of k or more samples with mutation for each given mutated site. The Poisson binomial calculations are made possible using the [poibin R pakage](https://cran.r-project.org/web/packages/poibin/poibin.pdf).

# Usage

python Main.py --M MutationFileListFile --C CovariateFileListFile --CC CombinedCovariateFile --LR LRModelName --P parallel --MF MergedMutationFilename --G grid --L logFilePath --RS regionSize

|Option | Description|
|--- | ---|
|MutationFileListFile | This should be a tab delimited file with patient id, mutation file location, and additional info (see below).|
|CovariateFileListFile | This should be a tab delimited file with covariate file name and file location.|
|CombinedCovariateFile | This should be a filename for the combined covariates. It can be generated as an intermediate but the name should be specified.|
|LRModelName | The name of the logistic regression model.|
|parallel | The number of jobs to run in parallel.|
|MergedMutationFilename | The name of the merged mutation file that is generated as an intermediate.|
|grid | 'T' to use grid engine. This option is not enabled yet.|
|logFilePath | The path to a log file (only used if grid option is 'T')|
|regionSize | Optional region size. Right now '1' is the only acceptable input.|

## Description of MutationFileListFile
This file should contain the following columns: pid,	MutationFile,	MutationWigFile,	MutationCovariateFile,	CoverageWigFile,	WGCovariateFile,	MutationCovariateSummaryFile,	ModelData

'pid' is the patient id. The mutation file is an input file the others are intermediate files generating during the application run.


# Steps

## Get Covariates for Mutations
These are base pair (AT or CG), replication timing, and coding/noncoding exon/intron annotations from GENCODE.

## Get Covariates for the Whole Genome
Same as for mutations but accross all bases with high coverage in the original WGS sequencing.
   
## Generate Sample Specific Probability Model Using Logistic Regression
Fit the logistic regression model to all data from all samples.
   
## Get Mutation Counts for Each Mutation Across Samples
Get the numbers of times a given genomic position is mutated across samples.
   
## Get Sample Specific Probabilities
Use the sample specific probability model from above to compute the site specific mutation probabilities for each sample.
   
## Compute Poisson Binomial Recurrence Probabilities
Using a vector of probabilities (one for each sample) and the poibin package compute the probabilities of seeing the observed number of mutations at each given site.



