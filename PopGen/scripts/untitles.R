library(RcppCNPy)

setwd("C:/Users/regina/Desktop/githubrepo/eco-genom/PopGen/results")

list.files()

#read in selection stats

s <- npyLoad("allRS_poly.selection.npy")
#rows are each locus, and 1 and 2 are the 2 eigan value groups with chisquare stats that we can convert to p-values

pval <- as.data.frame(1-pchisq(s,1))

head(pval)
#smaller pval shows selection

names(pval) = c("p_PC1", "p_PC2")
head(pval)

#associate those p-vals with our SNP meta data

p <- read.table("allRS_poly_mafs.sites", sep="\t", header = T, stringsAsFactors = T)

dim(p)
head(p)

p_filtered = p[which(p$kept_sites==1),]
dim(p_filtered)

cutoff=1e-3





## make manhattan plot
plot(-log10(pval$p_PC1),
     col=p_filtered$chromo,
     xlab="Position",
     ylab="-log10(p-value)",
     main="Selection outliers: pcANGSD e=1 (K2)")

#so our outliers are the upper log10, so these are the high outliers of pc1 showing BS hybridization
#log10 of 4 = 0.0001 > p 

#which contig has the really strong evidence of selection by  looking at p filtered
sel_contig <- p_filtered[which(pval==min(pval$p_PC1)), c("chromo", "position")]
sel_contig



outlier_contigs <- p_filtered[which(pval$p_PC1<cutoff),c("chromo", "position")]
head(outlier_contigs)
outlier_contigs

outlier_contigs <- outlier_contigs[which(outlier_contigs$position>0),]
dim(outlier_contigs)

outlier_contigs2 <- p_filtered[which(pval$p_PC2<cutoff),c("chromo", "position")]
dim(outlier_contigs2)
outlier_contigs2


write.table(unique(outlier_contigs$chromo),
            "allRS_poly_PC1_outlier_contigs.txt",
            sep="\t",
            quote=F,
            row.names=F,
            col.names = F)



write.table(unique(outlier_contigs2$chromo),
            "allRS_poly_PC2_outlier_contigs2.txt",
            sep="\t",
            quote=F,
            row.names=F,
            col.names = F)
