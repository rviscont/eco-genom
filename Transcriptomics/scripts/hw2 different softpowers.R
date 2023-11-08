#day 7 transcriptomics- keep it going

## Set your working directory
setwd('C:/Users/regina/Desktop/githubrepo/eco-genom/Transcriptomics/data')

# Load the package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

library(DESeq2)
library(ggplot2)

library(tidyverse)

install.packages("remotes")
remotes::install_github("kevinblighe/CorLevelPlot")

library(CorLevelPlot) 
library(gridExtra)

library(Rmisc) 
BiocManager::install("impute")

###################################################################################

#importing our data

################################################################################### 
# 1. Import the counts matrix and metadata and filter using DESeq2
#she wanted the filtered assembly so we could limit our input for wgcna
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
# conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
# head(conds)

sample_metadata = read.table(file = "Ahud_trait_data.txt",header=T, row.names = 1)

#a different colData from yesterday
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=sample_metadata, 
                              design= ~ 1)

dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# [1] 25260, that have at least 15 reads (a.k.a counts) in 75% of the samples

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds) 

###################################################################################

#outlier detection :)

###################################################################################

# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes
#wgcna has a gsg function
gsg <- goodSamplesGenes(t(countsTable))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)


# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(countsTable)), method = "average")
plot(htree) 
#DENDROGRAM tells you the outliers but we know that these aren't really from the gsg function
#but it does identify our weird groups as out groups


# pca - method 2

pca <- prcomp(t(countsTable))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

PCAtest <- ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
#we dont have to exclude anything here

###################################################################################

#3. Normalization ----------------------------------------------------------------------

################################################################################### 


colData <- row.names(sample_metadata)

# making the rownames and column names identical
all(rownames(colData) %in% colnames(countsTableRound)) # to see if all samples are present in both
all(rownames(colData) == colnames(countsTableRound))  # to see if all samples are in the same order



# perform variance stabilization
dds_norm <- vst(dds)
# dds_norm <- vst(normalized_counts)

# get normalized counts aka assay this and then transform it
norm.counts <- assay(dds_norm) %>% 
  t()
#################################################################################

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers

###################################################################################
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function; this step takes a couple minutes, signed means pick up or down regd, unsigned means absolute value 
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)
# based on this plot, choose a soft power to maximize R^2 (above 0.8) and minimize connectivity
# for these ahud data: 6-8; Higher R2 should yield more modules.


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 8 #SHOULD TRY RUNNING WITH 7 BASEDD ON RSQUARFED
temp_cor <- cor
cor <- WGCNA::cor # use the 'cor' function from the WGCNA package


# this step also takes a few minutes; ideally your maxBlockSize is larger than your number of genes to run the memory-intensive network construction all at once.
#Can minimize size of block size if you lack computational power, but its bigger than our # of genes so it'll use the whole data set
#to call something a module we set that it needs to have at least 30 genes (the standard)
#
bwnet8 <- blockwiseModules(norm.counts,
                          maxBlockSize = 26000,
                          minModuleSize = 30, 
                          reassignThreshold=0,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = F,
                          randomSeed = 1234,
                          verbose = 3)

# TOMtype (Topological Overlap Matrix type) parameter - unsigned - doesn't consider positive/negative co-expression
# signed - when you want to consider the direction of co-expression interaction, e.g., activating or inhibiting
# WGCNA often uses a dendrogram-based approach to identify modules. The choice of the 
# height cut in the dendrogram can determine the number of modules. Selecting a higher
# cut height results in fewer, larger modules, while a lower cut height leads to more, 
# smaller modules. so more than 25% difference b/n genes puts them in different modules

#bwnet has a lot of the important takeaways for visualization
#see 11 modules with gray being the random/non-caught genes
#could increase the soft-power threshhold to breakup the large groups
#theres more code to run in this tutorial we just didn't get there. 

#dendrogram shows you the clustering- again lots of turq and blue

cor <- temp_cor


bwnet8 <- readRDS("bwnet8.rds")
bwnet10 <- readRDS("bwnet10.rds")
bwnet14 <- readRDS("bwnet14.rds")
bwnet18 <- readRDS("bwnet18.rds")


###################################################################################

# 5. Module Eigengenes #8 ---------------------------------------------------------

###################################################################################

module_eigengenes8 <- bwnet8$MEs

head(module_eigengenes8)


# get number of genes for each module
table(bwnet8$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet8$dendrograms[[1]], cbind(bwnet8$unmergedColors, bwnet8$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# grey module = all genes that doesn't fall into other modules were assigned to the grey module
# with higher soft power, more genes fall into the grey module


###################################################################################

# 6A. # 8 Relate modules to traits --------------------------------------------------
# module trait associations

###################################################################################

traits <- sample_metadata[, c(5,8,11,14,17)]
#these are all the phenotype fitness and such data leaving behind the CIs
#gonna pull the pheno data and see if these modules match/correlate with any of these traits


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

#this is calc of pearsons correlation
module.trait.corr8 <- cor(module_eigengenes8, traits, use = 'p')
module.trait.corr.pvals8 <- corPvalueStudent(module.trait.corr8, nSamples)



# visualize module-trait association as a heatmap

heatmap.data8 <- merge(module_eigengenes8, traits, by = 'row.names')

head(heatmap.data8)

heatmap.data8 <- heatmap.data8 %>% 
  column_to_rownames(var = 'Row.names')


names(heatmap.data8)

CorLevelPlot(heatmap.data8,
             x = names(heatmap.data8)[16:20],
             y = names(heatmap.data8)[1:15],
             col = c("blue1", "skyblue", "white", "pink", "red"))

#see sign in the grey module with mean fitness so likely a new group could be teased apart
#dev is development time
######################################################################################

# 5. Module Eigengenes #10---------------------------------------------------------

###################################################################################

module_eigengenes10 <- bwnet10$MEs

head(module_eigengenes10)


# get number of genes for each module
table(bwnet10$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet10$dendrograms[[1]], cbind(bwnet10$unmergedColors, bwnet10$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# grey module = all genes that doesn't fall into other modules were assigned to the grey module
# with higher soft power, more genes fall into the grey module


###################################################################################

# 6A. Relate modules to traits --------------------------------------------------
# module trait associations

###################################################################################

traits <- sample_metadata[, c(5,8,11,14,17)]
#these are all the phenotype fitness and such data leaving behind the CIs
#gonna pull the pheno data and see if these modules match/correlate with any of these traits


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

#this is calc of pearsons correlation
module.trait.corr10 <- cor(module_eigengenes10, traits, use = 'p')
module.trait.corr.pvals10 <- corPvalueStudent(module.trait.corr10, nSamples)



# visualize module-trait association as a heatmap

heatmap.data10 <- merge(module_eigengenes10, traits, by = 'row.names')

head(heatmap.data10)

heatmap.data10 <- heatmap.data10 %>% 
  column_to_rownames(var = 'Row.names')


names(heatmap.data10)

CorLevelPlot(heatmap.data10,
             x = names(heatmap.data10)[14:18],
             y = names(heatmap.data10)[1:13],
             col = c("blue1", "skyblue", "white", "pink", "red"))

#see sign in the grey module with mean fitness so likely a new group could be teased apart
#dev is development time
######################################################################################
######################################################################################

# 5. Module Eigengenes #14---------------------------------------------------------

###################################################################################

module_eigengenes14 <- bwnet14$MEs

head(module_eigengenes14)


# get number of genes for each module
table(bwnet14$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet14$dendrograms[[1]], cbind(bwnet14$unmergedColors, bwnet14$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# grey module = all genes that doesn't fall into other modules were assigned to the grey module
# with higher soft power, more genes fall into the grey module


###################################################################################

# 6A. Relate modules to traits --------------------------------------------------
# module trait associations

###################################################################################

traits <- sample_metadata[, c(5,8,11,14,17)]
#these are all the phenotype fitness and such data leaving behind the CIs
#gonna pull the pheno data and see if these modules match/correlate with any of these traits


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

#this is calc of pearsons correlation
module.trait.corr14 <- cor(module_eigengenes14, traits, use = 'p')
module.trait.corr.pvals14 <- corPvalueStudent(module.trait.corr14, nSamples)



# visualize module-trait association as a heatmap

heatmap.data14 <- merge(module_eigengenes14, traits, by = 'row.names')

head(heatmap.data14)

heatmap.data14 <- heatmap.data14 %>% 
  column_to_rownames(var = 'Row.names')


names(heatmap.data14)

CorLevelPlot(heatmap.data14,
             x = names(heatmap.data14)[12:16],
             y = names(heatmap.data14)[1:11],
             col = c("blue1", "skyblue", "white", "pink", "red"))

#see sign in the grey module with mean fitness so likely a new group could be teased apart
#dev is development time
######################################################################################
######################################################################################

# 5. Module Eigengenes #18---------------------------------------------------------

###################################################################################

module_eigengenes18 <- bwnet18$MEs

head(module_eigengenes18)


# get number of genes for each module
table(bwnet18$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet18$dendrograms[[1]], cbind(bwnet18$unmergedColors, bwnet18$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# grey module = all genes that doesn't fall into other modules were assigned to the grey module
# with higher soft power, more genes fall into the grey module


###################################################################################

# 6A. Relate modules to traits --------------------------------------------------
# module trait associations

###################################################################################

traits <- sample_metadata[, c(5,8,11,14,17)]
#these are all the phenotype fitness and such data leaving behind the CIs
#gonna pull the pheno data and see if these modules match/correlate with any of these traits


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

#this is calc of pearsons correlation
module.trait.corr18 <- cor(module_eigengenes18, traits, use = 'p')
module.trait.corr.pvals18 <- corPvalueStudent(module.trait.corr18, nSamples)



# visualize module-trait association as a heatmap

heatmap.data18 <- merge(module_eigengenes18, traits, by = 'row.names')

head(heatmap.data18)

heatmap.data18 <- heatmap.data18 %>% 
  column_to_rownames(var = 'Row.names')


names(heatmap.data18)

CorLevelPlot(heatmap.data18,
             x = names(heatmap.data18)[8:12],
             y = names(heatmap.data18)[1:7],
             col = c("blue1", "skyblue", "white", "pink", "red"))

#see sign in the grey module with mean fitness so likely a new group could be teased apart
#dev is development time
######################################################################################
