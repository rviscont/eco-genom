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

#Entry 9: 2023-10-11
-f0 being egg to adult; F2 being 3 gens
-38 samples of F/R so 76 reads 

What questions/hypotheses can we test with this experimental design?
1. What are the effect of environment (temp and pH) on gene expression? Are there DEGs between all 4 treatment groups? What is unique to treatment? What is synergystic? 
2. What signatures of ( physiological ) adaptation and differences in gene expression of later generations exist? {compare F2, F4, and F11}
3. Does the gene expression in a given environment alter between generations? Within each environment and across environments? (interaction b/n enviro and generation)
10. Related to the above questions, what are the functional gene annotation catergories? 

Pipeline: 
  -clean raw seq with FastP - we will try to leave 'em long instead of trimming unknowns
  -Melissa already used Trinity to de novo assembly the ref transcriptome and asses it with BUSCO for completeness
  -use Salmon to map reads to the ref transcriptome but instead of hard mapping it'll use k-mers to map and quantify abundance
  -import into R for DESeq2 testing and visualization -> this si where we will spend a lot of our time in this unit
  
13 groups of treatment
everybody takes 1 group: I am doing Acidification (AH) F0 with 3 replicates

sign into our server; see our files with reps as the bucket they came from and the 1=F and 2=R
zcat AH_F0_Rep1_1.fq.gz | head -n 8
get to see Qs all between 20-40 so 99-99.99% confidence

cp /data/project_data/RNAseq/scripts/fastp_ahud.sh ~/myscripts/
get our script of interest
vim to edit
in the html file you can see the fastp version which is important to note within your hw/writeups
other important details: total reads before and after filtering, and how many reads passed doc'd in the good sheet

most important in the html files is the messiness in the base content ratios that is the stuff we didn't trim


#Entry 10: 2023-10-16 (Day 3)
TODAY: check our google sheet, ass our de novo assembly quality, mapping and quantifying abundance using Salmon, 1 person will pool our samples into 1 table which we will then use to input DESeq2 to analyze gene expression

-google sheet shows average # of reads, average # kept, and % retained - A GOOD VALUE TO REPORT
-keep in mind the data you generate depends on the samples you put in - if you are looking for low expression genes you might need more samples/more stratification specific to what you are looking for (tissue type, age, etc.)
-more data means more contigs for trinity de novo assembly 

USE TRINITY TO ASSESS QUALITY
/data/popgen/trinityrnaseq-v2.13.2/util/TrinityStats.pl  /data/project_data/RNAseq/assembly/ahud_Trinity.fasta
-uses pearl language to open the .fasta
-reminder N50 is the string of all contigs and the middle length is that number - longer than the avg and median contigs cuz lots of little seqs get tacked onto the end of the string of all contigs so bring down avg/median
-if i was gonna test just 1 rep per treatment group, aka 1/3 of the data, most of these measurements would improve- cuz we gave it a lot of data, and its counting isoforms/splice variants/allele variants as separate genes, we are gonna filter this down :) 

BUSCO
-looking for single copy orthologs
busco -m transcriptome -i /data/project_data/RNAseq/assembly/ahud_Trinity.fasta -o ~/myresults/BUSCO/BUSCOarthropoda -l arthropoda_odb10
-assessing transcriptome mode, input is assembly, output busco directory, then compare to arthropods
-creates summary shown in the tutorial day 3
-complete buscos=96.9%/1013 orthologs
-not surprising there were a lot of duplicates with our giant data set
-to make ref transcriptome, take the 100 single copies, and deleting the duplicates: we are going to keep the duplicates to help assess potential variation

MAPPING REF TRANSCRIPTOME by SALMON
-Salmon is part of trinity, using pearl language
-prep the data, 1 person has to do so for whole class

/data/popgen/trinityrnaseq-v2.13.2/util/align_and_estimate_abundance.pl --transcripts /data/project_data/RNAseq/assembly/ahud_Trinity.fasta \
  --est_method salmon \
  --trinity_mode \
  --prep_reference
  
- makes 2 index files

Mapping the Ref using Salmon
-ahud_samples.txt > ahud_AH_F0.txt

cd /data/project_data/RNAseq/mapping

/data/popgen/trinityrnaseq-v2.13.2/util/align_and_estimate_abundance.pl --transcripts /data/project_data/RNAseq/assembly/ahud_Trinity.fasta \
  --seqType fq \
  --samples_file ~/mydata/ahud_AH_F0.txt \
  --est_method salmon \
  --output_dir /data/project_data/RNAseq/mapping \
  --thread_count 1 \
  --trinity_mode
  
at the end we can do this: to get mapping rate
grep -r --include \*.log -e 'Mapping rate' > ~/myresults/ahud_mapping_rate

whatever we get is the percent of reads that mapped to the 3 thousand something 'genes'


#Entry 11: 2023-10-18 (Day 4)
-check our mapping
-can be affected by read quality (but we already cleaned ours), 
-why might that one sample of 68% not be mapping well? well could be not a copepod, might blast it, something to think about....

-in mapping folder you can go into each treatment rep - 
-in quant.sf TPM is how many transcripts per million and num_reads being unscaled
-then we changed the label of each file to its base directory (treatment/sample)
-creating script for today:counts data
-grabbed *.matrix files that were just generated which are the counts data
-38 samples, ~30thou isoforms
-rounding cuz DESeq doesn't like decimals
-conds is our conditions
-

-exploring the count data, we end up with a mean of ~20million reads?
-szie factors is caling
-dispersion estimates is how dispersed a given gene; reads ranging from 10-1000 means big dispersion vs always highly expressed small dispersion;
-DESeq fit our model and did our testing for us
-

### Entry 12: 2023-10-23 (Day 5).
-today we take a step back and make sure that the decisions we made were good ones
-used cdkit and transdecode to get rid of allele variants (95% similarity) and length variants
-see our new stats on mapping, and busco scores in the tutorial
-am amazing!!! weightedd deg venn diagram- gene specific comparison across treatment or generation
-also heat maps

#Entry 13: 2023-10-25 (Day 6)
-today we use wgcna aka gene network analysis
-weighted gene co-expression network analysis
-it will look for clusters called modules between transcripts across samples
-based on matrix algebra and eigengenes basically describing the pattern across samples
-processing our salmon table to play in wgcna then determining the most appropriate amount of modules (power analysis and we choose the Rsquared we desire dedpendng on how tight we want each group) and finally visualizing the data 
-wgcna finds modules independent of the samples so they cluster across treatments (agnostic to treatments)
-then we gived it the trait metadata aka the physiological data and ask if theres a correlation between our model and phenotypes
-

### Entry 14: 2023-10-30.
-more WCGNA analysis day 7 wgcna script

-transferred file from server after Melissa
  "# What I did to save it:
  saveRDS(bwnet, file = "bwnet.rds")
  # To load the object
  bwnet <- readRDS("bwnet.rds")"
  
-re-ran code doing counts and normalization, but not outlier detection or network construction cuz Melissa provided the bwnet file
-eigengene - negative = low expression of the genes in that module in that sample
-increase power threshhold (from 6 -> 10 maybe) you'll get a higher Rsquared so number of module groups will increase -- good if we wanted higher res of these groups play with power threshold
-for figure, the merged is collapsed modules from the unmerged groups
-line interaction plot shows treatment groups across generations - this is the pattern of the eigengene expression - so these genes were grouped in the yellow module because they are kinda doing the same thing - again see difference b/n OWA and OW vs OA and AM -
-tophub genes show the same thing as this prev plot and more - lots of similarity in gen 4
-heat map is the 91 genes making up the purple module
-not a whole lot of clustering by gen or by treatment
-changed black to white on heat maps and saw more consistant clustering
-2 classes of GO - 1) a significant cutoff (FET) and compare sig vs non-sig usually by Phishers or chi-square so it'll go thru each category making a 2x2 of sig and non sig x vs gene and non gene on y
-if theres no enrichment (selection) you'll see the same ratios of the of total genes in the table
-if there IS enrichment, we look at observed vs expected and that tells us there is significance 
-2) a distribution approach / gene set enrichment analysis using rand-based MVU or KS test- take the -log(pval) and graph, where very positive is sig, per each gene category they search and find that the the distribution is different from the genome as a whole

### Entry 15: 2023-11-01.
-day 8 DESeq2_to_GOMWU.R
-hw review: pick a question, and see how techincal decisions affect bio results
-5 options:
1 filtration- we started with the whole assembly, but then filterered it; but maybe we should just use cd-hit to cluster by similarity instead of adding the transdecoder program
2 within DEseq we filtered by read depth, but we could play with how depth coverage might affect our transcriptome
3 within wgcna change the softpower threshhold to see the module formation, correlation, traits, with the bonus of adding the GO enrichment which is used to validate what you find with the modules
4 bio questions: not technical-  gene expression b/n AM and OWA acorss generations; also GO enrichment
5 bio questions: not technical-differences in gene expression between F0 and F4; also GO enrichment- can use this to answer if OWA suggests synergy, additon or antagonization of the gene patterns from OW and OA
6 maybe another quesiton.....???
-same criteria as before, but up to 4 figures/tables also as both as word doc and PDF
-take comments from steve into account
-due Tuesday end of day

-we get to say min/max for GO analyses
-10 genes to have statistical power
-max at 500 genes so we can stay at the level of meaningful categories
-we can do enrichment on any quantitative measurement, or on wcgna results (categorical) - usually a distribution
-be aware the doc we are using are treating isoforms like totally different genes
-

### Entry 16: 2023-11-06.

-Csenge's lesson on structural variation
-echo is the command line
-cd /netfiles/ecogen/structural_variation/
-im chromosome 5 -  NW_022145598.1
-filter_chromosome.sh - we make a vcf thats convelrted back to bcf to be indexed to speed up searching.
-vcf tells you how each genetic variant is different from the ref per 1 chromosome across all 140 individuals
-head tells you how the file was made, and then about the ref genome, then it goes to the list 
-we use more -line# to see further down, hit enter to continue to navigate
-bottom of file tells you what all the codes mean.
-format is GT=genotype and PS=phred score
-not just biallelic, and for each gene you see ref=0 and the options of what it is
-period means not known, 
-liklihoods for each genotype listed ":X,Y,Z" with number of liklihoods depending on the # of different allele combos - and biggest number is less likely 

