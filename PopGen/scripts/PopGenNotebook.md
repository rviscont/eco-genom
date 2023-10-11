# Population Genomics Lab Notebook

## Author: Regina Visconti

### Affiliation: University of Vermont: Department of Plant Biology; Department of Cellular, Molecular, and Biomedical Sciences

### E-mail contact: rviscont\@uvm.edu

### Start Date: 9/11/23

### End Date:

### Project Descriptions: This notebook is to document my bioinformatics workflow for Population Genomics Fall 2023.

# Table of Contents:

-   [Entry 1: 2023-09-11](#id-section1)
-   [Entry 2: 2023-09-13](#id-section2)
-   [Entry 3: 2023-09-18](#id-section3)

|     |
|-----|
|     |
|     |

### Entry 1: 2023-09-11.

-   Reviewed red spruce study system exome capture data
-   Discussed fastq files (DNA seq +q-score)
-   used FASTQC to analyze sequencing runs

### Entry 2: 2023-09-13.

-   FASTQC revealed high quality data for read length
-   eliminated any wonky reads (low Q scores)
-   developed read trimming code by fastp
-   made html files to compare pre- and post- trim
-   set up read mapping of trimmed/clean reads using 'bwa'

### Entry 3: 2023-09-18.

- Introduce lab notebooks
- Visualize sequence alignment files (*.sam)
- Process the mapping file *sam to binary (*.bam), sort, and remove duplicate reads from PCR
- Calculate mapping statistics to assess quality of the result


- 
### Entry 6: 2023-09-27.
-looking at FST to ask how much gene flow is between black spruce and our pops
-Fst is variance in allele frequencies =variance(minor allele)/avg allele freq
-gonna use ANGSD site freq spectrum SFS to compare two SFS
-Fst is usually 0-0.05 for 1 plant species
-human Fst is about the same
-Fst urchin ~0.01
-large verts like whales or migratory birds Fst goes up with range ~0.1
-b/n sister taxa like blackspruce and our pop itll be ~0.1-0.2



-looking at our google doc theta pi increases going noth
-same for theta(w) segregating sites aka how many of the bps were polymorphic
-seems to be affected by isolation in the south
-Tajimas D 0 is equilibrium, population isn't growing or shrinking
-if D distribution shifts positive, it means loss of rare alleles aka a bottleneck 
-if D shifts negative, it means an abundance of rare alleles (Ne * mutation) so usually means a growing pop
-so for ours the southern pops are more positive as the pops are small due to fragmentation and isolation
-can estimate Ne from theta 
-theta =  4 * Ne * mutation per bp per generation time [30 years]
-4 cuz theyre diploid (2 copies) and 2 parents
-so with theta w and pi solve fro Ne to get a range
-Ne = theta/4* mutation(30x10^-9)



-wrote a script to make a PCA plot and admixture analysis
-so we will get our clusters and estimations of allele frequencies for each individual in each pop
-genotype liklihood is input and take a long time to calc
-the script is given on day 6 on the popgen tutorials

#Entry 7: 2023-10-02

Reviewing Plot and Admix Plots
-Eigenvalues plot: a scree plot derived from the PCA
  -PCA establishes a linar model between the SNPs to best explain the variation in the allele frequencies = PC1
  - then our model slices through the cloud of NSPs orthoganally over and over to get our following PCs
  -so the eigenvalues show how much variability is accounted for by each PC
  
-the pretty PCA plot
  -see 2021 is pretty isolated from the rest of the pops by PC2 -> PA pop where there are no high mountains, kinda boggy, boreal enclave 
  -2020 and 2019are sorted along both axes -> these are the southern
  appalacia pops
  -along the PC1 axis, see some sorting of 2100 and 2101, so maybe its gene flow south to north
  
-admixture plot
  -each thin bar is an individual
  -y-axis "admix proportions" which means that if we split our pops into 2, 1 is blue 1 is orange, can see mixing of the pops F1s be ~50/50
  -mixing can indicate backcrossing, and we kinda see that 2022/24/27
  -overall see increasing freq of blue going north, probs driven by hybridization with Northern Black Spruce
  
-Compare to our Fst values
  -its for each pop vs the black spruce sample (which were far away from our red spruce pops)
  -see decreasing value as you move more north
  -lower Fst shows more hybridization so this is consistent with our data
  -2021 comes out again as pretty separate from everyone else since it has a lower Fst so kinda weird it is hybridizing even though it is so isolated
  -also confirms our diversity values too

  -admix plot and PCA are from the same data so we can see similar trends -> admix kinda shows PC1 so maybe PC1 is hybridizaiton (aka a type of gene flow)
  
  
So now we can think about selection
-we have to specifically test for it even though we've hinted at it
-section 2 of daily notes pop gen 7.1
-semi-normal distribution of sensitive genetics like maybe Fst by freq
-so our outliers are the tails of the distribution so same data, different plot, trying to establish whats drift/neutral processes and with natural selection drives us to extreme Fsts
-high values=local selection and adaptation (divergence between pops)
-low Fst = allele freqs are shared among pops so if you have an allele that's beneficial everywhere, so balancing selection OR universally favored that aren't fixed yet.
-USE PCangsd.SH to test for these outliers

-did a lot of R and more; 
-going from outlier loci, to which contigs theyre on, to which genes they are!
-now we can use a database to annotate them (what they do) and if they are enriched for certain functions

#Entry 8: 2023-10-04

-gene-environment association analysis 1 locus at a time
-need genotype data (our SNPs) and environment data 
-we will use covariates to make sure we are looking at the right things usually a measure of co-pop structure like PC axes or admix coefficients 
-covariates are things that might change between genes/environment that aren't causal 

-pipeline: starting in R
  1) analysis of outlineon outlier loci and export the list specifically related to PC1
  2) genetic PCA and export PC1 and PC2 as our co-variates
  3)get bioclim(environment) variables and export them to test
  4)transfer files to server to run GEA
  
#Entry 9: 2023-10-09

FIRST DAY TRANSCRIPTOMICS

-geno causes pheno in context of the environment
-so pop geno is looking at SNPs and environment aka gneo and environment and can be intergrated with pheno
-1 pheno is the whole transcriptome; all genes being expressed
-things that affect pheno: enviro conditions (treatment, precipitation/temp/bioclim variables[which are long term] or local conditions for field sample transcriptomics, acclimation period or common garden where we try to hold enviro constant), time (of year, of day) /devo stage, genetic background/population/sex, and tissue type sample
-transcriptomics are a pheno snap shot and all these variables that affect pheno need to be included as part of experimental design [depending on your ?/hypothesis]; you are making decisions that affect your outcomes the whole way
-pt of the experimental design is extracting RNA; then what kind of library do you need to prep/seq? then make cDNA
      +you can do polyA tail prep which selects for the mRNA with the polyA tail meant to direct them to transcription/translation/expression;
      +miRNA likely to interact with more than 1 mRNA strand and so can pt out some big regulators
-sequence!
      +Illumina (everyone's favorite for the short mRNA transcripts unless looking for splice variants)
      +Ion Torrent
      +Oxford Nanopore (ONT)
      +PacBio (ONT and PB are good for longer splice variant tscs)
-now you have data dn can do ANALYSIS AND VISUALIZATION
      +FASTQ files (fasta + quality) used to visualize and clean data to get clean_fq files [we use fastp]
      + if not a MO, may use clean_fq to do a de novo assembly so you can map: itierative process, eval,[using Trinity] - biggest issue w/ this step is differentiating b/n splice variants/paralogs/orthologs so they just call them isoforms
      +map these reads to published/de novo assembly to get read counts [using parental salmon]
      +then normalization [DESeq2]
      +identify the DEGs [DESeq2]
      +now you can do GO, or Weighted Gene Correlation Network Analysis on the outputs [we'll use GOMWU or TOPGO; WGCNA]
-we aren't looking at nts unless trying to assemble; we are looking for # of reads
-