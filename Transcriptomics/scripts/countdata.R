#counts data for RNA DESeq2

#setwd('C:/Users/regina/Desktop/githubrepo/eco-genom/Transcriptomics/data')

#install.packages("BiocManager")

#BiocManager::install("DESeq2")

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)  ### First: BiocManager::install("vsn") AND BiocManager::install("hexbin")

###################################################################################

#importing our data

###################################################################################

# Import the counts matrix
countsTable <- read.table("salmon.isoform.counts.matrix", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample discription table
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)


###################################################################################

#explore count data and distributions

###################################################################################

# Let's see how many reads we have from each sample
colSums(countsTableRound)
mean(colSums(countsTableRound))

barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound),cex.names=0.5, las=3,ylim=c(0,21000000))
abline(h=mean(colSums(countsTableRound)), col="blue", lwd=2)

# the average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) # [1] 2245.401
median(rowSums(countsTableRound)) # [1] 117

apply(countsTableRound,2,mean) # 2 in the apply function does the action across columns
apply(countsTableRound,1,mean) # 1 in the apply function does the action across rows
hist(apply(countsTableRound,1,mean),xlim=c(0,1000), ylim=c(0,300000),breaks=10000)
#distribution shows a few genes with high counts, and a lot of genes with not a lot of counts?

###################################################################################

#start using DESeq2!

###################################################################################

#### Create a DESeq object and define the experimental design here with the tilda according to DESeq2 tutorial conventions

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ generation + treatment)

dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ-- we picked 30 so that we have fewer isoforms to analyze, this is one of our many decisions

dds <- dds[rowSums(counts(dds) >= 30) >= 28,]
nrow(dds) # 41348; COULD ALSO REPlot at this time.


# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)
#ask your question cuz DESeq can't really test interaction across so many treatments/samples, so know what you want to know. 