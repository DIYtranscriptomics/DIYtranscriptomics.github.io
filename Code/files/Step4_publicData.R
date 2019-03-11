# Introduction to this script -----------
# the goal of this script is to give you access to massive amounts of gene expression data without the need to download and analyze individual .fastq files
# this script allows you to access RNAseq data from GEO and SRA using the "all RNA-seq and ChIP-seq sample and signature search" (ARCHS4) project (https://amp.pharm.mssm.edu/archs4)
# this script also provides programmatic access to the "Library of Integrated Network-Based Cellular Signatures" (LINCS)  

# Load packages ------
library(tidyverse)
library(reshape2)
library(rhdf5)
library(edgeR)

# load ARCHS4 database -----
# you should have already downloaded the most recent versions of mouse and human RNAseq data from ARCHS4 in hdf5 format
# begin by setting the path to your archs4 data files
archs4.human <- "~/Dropbox/publicData/human_matrix.h5"
archs4.mouse <- "~/Dropbox/publicData/mouse_matrix.h5"
# use the h5 list (h5ls) function from the rhdf5 package to look at the contents of these databases
# but first, remember our Kallisto outputs were hdf5 files, so let's look at one of these first.
h5ls("../readMapping/uninf_rep1/abundance.h5")
h5ls(archs4.human)
h5ls(archs4.mouse)

# 133,776 samples from human
all.samples.human <- h5read(archs4.human, name="meta/Sample_geo_accession")

# 170,010 samples from mouse
all.samples.mouse <- h5read(archs4.mouse, name="meta/Sample_geo_accession")

# query ARCHS4 database ----
# choose your samples based on GEO or SRA ID
mySamples <- c("GSM2310941", # WT_unstim_rep1
               "GSM2310942", # WT_unstim_rep2
               "GSM2310943", # Ripk3_unstim_rep1
               "GSM2310944", # Ripk3_unstim_rep2
               "GSM2310945", # Ripk3Casp8_unstim_rep1
               "GSM2310946", # Ripk3Casp8_unstim_rep2
               "GSM2310947", # WT_LPS.6hr_rep1
               "GSM2310948", # WT_LPS.6hr_rep2
               "GSM2310949", # Ripk3_LPS.6hr_rep1
               "GSM2310950", # Ripk3_LPS.6hr_rep2
               "GSM2310951", # Ripk3Casp8_LPS.6hr_rep1
               "GSM2310952") # Ripk3Casp8_LPS.6hr_rep2

# Identify columns to be extracted from ARCHS4 database
my.sample.locations <- which(all.samples.mouse %in% mySamples)
# extract gene symbols from the metadata
genes <- h5read(archs4.mouse, "meta/genes")

# Extract expression data from ARCHS4 ----
expression <- h5read(archs4.mouse, "data/expression", index=list(1:length(genes), my.sample.locations))
H5close()
rownames(expression) <- genes
colnames(expression) <- all.samples.mouse[my.sample.locations]
colSums(expression) #this shows the sequencing depth for each of the samples you've extracted
archs4.dgelist <- DGEList(expression)
archs4.cpm <- cpm(archs4.dgelist)

# Filter and normalize extracted data ----
table(rowSums(archs4.dgelist$counts==0)==9)
keepers <- rowSums(archs4.cpm>1)>=2
archs4.dgelist.filtered <- archs4.dgelist[keepers,]
dim(archs4.dgelist.filtered)
archs4.dgelist.filtered.norm <- calcNormFactors(archs4.dgelist.filtered, method = "TMM")

archs4.filtered.norm.log2.cpm <- cpm(archs4.dgelist.filtered, log=TRUE)
colnames(archs4.filtered.norm.log2.cpm) <- names

# Extract sample metadata from ARCHS4 to create a study design file ----
# extract the sample source
Sample_source_name_ch1 <- h5read(archs4.mouse, "meta/Sample_source_name_ch1")
# extract sample title
Sample_title <- h5read(archs4.mouse, name="meta/Sample_title")
# extract sample characteristics
Sample_characteristics<- h5read(archs4.mouse, name="meta/Sample_characteristics_ch1")

# let's try putting this all together in a study design file
studyDesign <- tibble(Sample_title = Sample_title[my.sample.locations], 
                   Sample_source = Sample_source_name_ch1[my.sample.locations],
                   Sample_characteristics = Sample_characteristics[my.sample.locations])

#based on what we extracted from ARCHS4 above, lets customize and clean-up this study design file
studyDesign <- tibble(Sample_title = Sample_title[my.sample.locations], 
                   genotype = c("WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8", "WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8"),
                   treatment = c("unstim", "unstim", "unstim", "unstim", "unstim", "unstim", "LPS", "LPS", "LPS", "LPS", "LPS", "LPS"))

# Pricipal component analysis (PCA) -------------
pca.res <- prcomp(t(archs4.filtered.norm.log2.cpm), scale.=F, retx=T)
#look at pca.res in environment
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
x <- pca.res$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
pca.res$x #$x shows you how much each sample influenced each PC (called 'loadings')
#note that these loadings have a magnitude and a direction (this is the basis for making a PCA plot)
pc.var<-pca.res$sdev^2 #sdev^2 gives you the eigenvalues
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per

# Visualize your PCA result ------------------
#lets first plot any two PCs aslgainst each other
#We know how much each sample contributes to each PC (loadings), so let's plot
pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df, aes(x=PC1, y=PC2, color=studyDesign$treatment)) +
  geom_point(size=4) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_ipsum_rc()

# now try painting other variables from your study design file onto this PCA.
# can you determine the relationship between PC2 and your metadata?

# now create a small multiple PCA plot
melted <- cbind(factor(studyDesign$treatment), melt(pca.res$x[,1:4]))
head(melted) #look at your 'melted' data
colnames(melted) <- c('group', 'treatment', 'PC', 'loadings')
ggplot(melted) +
  geom_bar(aes(x=treatment, y=loadings, fill=group), stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_flip() +
  theme_ipsum_rc()
