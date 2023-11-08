#day 6 transcriptomics

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
 
 ggplot(pca.dat, aes(PC1, PC2)) +
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
 
 # 4. Network Construction  ---------------------------------------------------
 # Choose a set of soft-thresholding powers
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
 
 soft_power <- 6 #SHOULD TRY RUNNING WITH 7 BASEDD ON RSQUARFED
 temp_cor <- cor
 cor <- WGCNA::cor # use the 'cor' function from the WGCNA package
 
 
 # this step also takes a few minutes; ideally your maxBlockSize is larger than your number of genes to run the memory-intensive network construction all at once.
 #Can minimize size of block size if you lack computational power, but its bigger than our # of genes so it'll use the whole data set
 #to call something a module we set that it needs to have at least 30 genes (the standard)
 #
 bwnet <- blockwiseModules(norm.counts,
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