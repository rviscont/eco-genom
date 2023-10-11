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

outliers_PC1_e2 <- p_filtered[which(pval$p_PC1<cutoff),c("chromo","position")]

# how many outlier loci < the cutoff?
dim(outliers_PC1_e2)[1]

write.table(outliers_PC1_e2,
            "allRS_poly_outliers_PC1_e2.txt", 
            sep=":",
            quote=F,
            row.names=F,
            col.names=F)


COV <- as.matrix(read.table("allRS_poly.cov"))

PCA <- eigen(COV)

data=as.data.frame(PCA$vectors)
data=data[,c(1:2)] # the second number here is the number of PC axes you want to keep

head(data)

write.table(data,
            "allRS_poly_genPC1_2.txt",
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)


#LAST POPULATION GENOMICS CODING SESH
#genotype environemnt association of selection outlier loci

#luckily we will be redoing pcANGSD scan for genetic PCs and outlier loci
#AND MORE



library(raster)#good for clim data
library(FactoMineR) #these factos are good for PCA in R
library(factoextra)
library(corrplot)#yay for pretty plots

bio <- getData("worldclim",var="bio",res=10)

coords <- read.csv("https://www.uvm.edu/~kellrlab/forClass/colebrookSampleMetaData.csv", header=T)
head(coords)
#The chunk below refers to your bamlist file that you transferred during last week's PCA/admixture analysis.  It should be the same one you want to use here -- if your sample list for analysis changes in the future, you'll need a different bamlist!
#theres 900 samples in the meta data, we just want the 95 we seqd based on the tree ID

names <- read.table("allRS_bam.list")
names <- unlist(strsplit(basename(as.character(names[,1])), split = ".sorted.rmdup.bam"))
split = strsplit(names, "_")
pops <- data.frame(names[1:95], do.call(rbind, split[1:95]))
names(pops) = c("Ind", "Pop", "Row", "Col")

#merged our cleaned up pops with metadata coords
angsd_coords <- merge(pops, coords, by.x="Ind", by.y="Tree")

points <- SpatialPoints(angsd_coords[c("Longitude","Latitude")])

clim <- extract(bio,points)

angsd_coords_clim <- cbind.data.frame(angsd_coords,clim)
str(angsd_coords_clim)



# Make the climate PCA: 15-33 cuz these are the bioclim variables in our file
clim_PCA = PCA(angsd_coords_clim[,15:33], graph=T)
#read it as arrows show correlation
#bio1 corrrelated to PC1 and 2 - annual mean temp so pops in that area of the graph have high temp and opposite that arrow the pops have low mean annual temp


# Get a screeplot of cliamte PCA eigenvalues - first PC1 explains more than 50% of the variability
fviz_eig(clim_PCA)

# What is the climate PCA space our red spruce pops occupy? easiest to interpret so far
fviz_pca_biplot(clim_PCA, 
                geom.ind="point",
                col.ind = angsd_coords_clim$Latitude, 
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                title="Climate PCA (Bioclim)",
                legend.title="Latitude")

# Which variables show the strongest correlation on the first 2 (COUNT EM 2) climate PC axes?
dimdesc(clim_PCA)[1:2]

#so we know the most important variable for each axes and we want to export that data
# Replace "XX" with your bio variable most significant on climate PC1:
write.table(scale(angsd_coords_clim["bio10"]),
            "allRS_bio10.txt",
            sep="\t",
            quote=F,
            row.names = F,
            col.names=F)


# Replace "YY" with your bio variable most significant on climate PC2:  
write.table(scale(angsd_coords_clim["bio12"]),
            "allRS_bio12.txt",
            sep="\t",
            quote=F,
            row.names = F,
            col.names=F)
