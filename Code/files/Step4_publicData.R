# Introduction to this script -----------
# the goal of this script is to give you access to massive amounts of gene expression data without the need to download and analyze individual .fastq files
# this script allows you to access RNAseq data from GEO and SRA using the "All RNA-seq and CHIP-Seq Sample and Signature Search" (ARCHS4) project (https://amp.pharm.mssm.edu/archs4)

# Load packages ------
# nothing new here...you should already have all these packages in your R package library
library(tidyverse)
library(rhdf5)
library(edgeR)

# load ARCHS4 database -----
# you should have already downloaded the most recent versions of mouse and human RNAseq data from ARCHS4 in hdf5 format
# begin by creating file paths that point to the hdf5 archs4 files
# because of the size of these files, feel free to skip the human data and just work with mouse
archs4.human <- "archs4_gene_human_v2.1.2.h5"
archs4.mouse <- "archs4_gene_mouse_v2.1.2.h5"
# use the h5 list (h5ls) function from the rhdf5 package to look at the contents of these databases
h5ls(archs4.human)
h5ls(archs4.mouse)

# 620,825 samples from human
all.samples.human <- h5read(archs4.human, name="meta/samples/geo_accession")
dim(all.samples.human)

# 717,966 samples from mouse!
all.samples.mouse <- h5read(archs4.mouse, name="meta/samples/geo_accession")
dim(all.samples.mouse)

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
my.sample.locations <- which(all.samples.mouse %in% mySamples) # first time you've seen the %in% operator.
# extract gene symbols from the metadata
genes <- h5read(archs4.mouse, "meta/genes/gene_symbol")

# Extract expression data from ARCHS4 ----
expression <- h5read(archs4.mouse, "data/expression",
                     index=list(my.sample.locations, 1:length(genes)))
# transpose to get genes as rows and samples as columns
expression <- t(expression)

rownames(expression) <- genes
colnames(expression) <- all.samples.mouse[my.sample.locations]
colSums(expression) #this shows the sequencing depth for each of the samples you've extracted
archs4.dgelist <- DGEList(expression)
archs4.cpm <- cpm(archs4.dgelist)
colSums(archs4.cpm)

# Filter and normalize the extracted data ----
table(rowSums(archs4.dgelist$counts==0)==12)
keepers <- rowSums(archs4.cpm>1)>=2
archs4.dgelist.filtered <- archs4.dgelist[keepers,]
dim(archs4.dgelist.filtered)
archs4.dgelist.filtered.norm <- calcNormFactors(archs4.dgelist.filtered, method = "TMM")

archs4.filtered.norm.log2.cpm <- cpm(archs4.dgelist.filtered.norm, log=TRUE)

# Extract sample metadata from ARCHS4 to create a study design file ----
# extract the sample source
sample_source_name <- h5read(archs4.mouse, "meta/samples/source_name_ch1")
# extract sample title
sample_title <- h5read(archs4.mouse, name="meta/samples/title")
# extract sample characteristics
sample_characteristics<- h5read(archs4.mouse, name="meta/samples/characteristics_ch1")

# let's try putting this all together in a study design file
studyDesign <- tibble(Sample_title = sample_title[my.sample.locations],
                      Sample_source = sample_source_name[my.sample.locations],
                      Sample_characteristics = sample_characteristics[my.sample.locations])

#based on what we extracted from ARCHS4 above, lets customize and clean-up this study design file
studyDesign <- tibble(Sample_title = sample_title[my.sample.locations],
                      genotype = c("WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8", "WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8"),
                      treatment = c("unstim", "unstim", "unstim", "unstim", "unstim", "unstim", "LPS", "LPS", "LPS", "LPS", "LPS", "LPS"))

#capture experimental variables as factors from this study design
genotype <- factor(studyDesign$genotype)
treatment <- factor(studyDesign$treatment)
sampleName <- studyDesign$Sample_title

# Principal component analysis (PCA) -------------
pca.res <- prcomp(t(archs4.filtered.norm.log2.cpm), scale.=F, retx=T)
#look at pca.res in environment
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
pca.res$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
pca.res$x #$x shows you how much each sample influenced each PC (called 'loadings')
#note that these loadings have a magnitude and a direction (this is the basis for making a PCA plot)
pc.var<-pca.res$sdev^2 #sdev^2 gives you the eigenvalues
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per

# Visualize your PCA result ------------------
#lets first plot any two PCs against each other
#We know how much each sample contributes to each PC (loadings), so let's plot
pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2) +
  geom_point(size=4) +
  # geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) +
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

# now try painting other variables from your study design file onto this PCA.
# can you determine the relationship between PC2 and your metadata?
# can we map one variable to point color and another to point shape?

# now create a small multiple PCA plot
pca.res.df <- pca.res$x[,1:4] %>% # note that this is the first time you've seen the 'pipe' operator from the magrittr package
  as_tibble() %>%
  add_column(treatment) %>%
  add_column(genotype) %>%
  add_column(sampleName)


pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca.pivot) +
  aes(x=sampleName, y=loadings, fill=treatment) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes. Try doing this with the genotype variable you created above.
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()



# the essentials ----
library(tidyverse)
library(rhdf5)
library(edgeR)

archs4.mouse <- "mouse_matrix_v10.h5" # if you placed the hdf5 file in your working directory, just use "human_matrix_v8.h5" as the path
all.samples.mouse <- h5read(archs4.mouse, name="meta/samples/geo_accession")
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

my.sample.locations <- which(all.samples.mouse %in% mySamples)
genes <- h5read(archs4.mouse, "meta/genes/gene_symbol")
expression <- h5read(archs4.mouse, "data/expression",
                     index=list(my.sample.locations, 1:length(genes)))

expression <- t(expression)
rownames(expression) <- genes
colnames(expression) <- all.samples.mouse[my.sample.locations]
archs4.dgelist <- DGEList(expression)
archs4.cpm <- cpm(archs4.dgelist)

keepers <- rowSums(archs4.cpm>1)>=2
archs4.dgelist.filtered <- archs4.dgelist[keepers,]
archs4.dgelist.filtered.norm <- calcNormFactors(archs4.dgelist.filtered, method = "TMM")
archs4.filtered.norm.log2.cpm <- cpm(archs4.dgelist.filtered.norm, log=TRUE)

sample_source_name <- h5read(archs4.mouse, "meta/samples/source_name_ch1")
sample_title <- h5read(archs4.mouse, name="meta/samples/title")
sample_characteristics<- h5read(archs4.mouse, name="meta/samples/characteristics_ch1")

studyDesign <- tibble(Sample_title = sample_title[my.sample.locations],
                      Sample_source = sample_source_name[my.sample.locations],
                      Sample_characteristics = sample_characteristics[my.sample.locations])

studyDesign <- tibble(Sample_title = sample_title[my.sample.locations],
                      genotype = c("WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8", "WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8"),
                      treatment = c("unstim", "unstim", "unstim", "unstim", "unstim", "unstim", "LPS", "LPS", "LPS", "LPS", "LPS", "LPS"))

pca.res <- prcomp(t(archs4.filtered.norm.log2.cpm), scale.=F, retx=T)
pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color=treatment, shape=genotype) +
  geom_point(size=4) +
  # geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) +
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

