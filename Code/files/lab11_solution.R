# Load packages ----
library(tidyverse)
library(rhdf5)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots

# get data from ARCHS4 -----
# load the sample IDs that you were given for today's lab
mySamples <- read_tsv("skin_ids.tsv")
mySamples <- mySamples$ids
# begin by creating file paths that point to the hdf5 archs4 files
archs4.human <- "human_gene_v2.3.h5"
all.samples.human <- h5read(archs4.human, name="meta/samples/geo_accession")
# Identify columns to be extracted from ARCHS4 database
my.sample.locations <- which(all.samples.human %in% mySamples) # first time you've seen the %in% operator.
# extract gene symbols from the metadata
genes <- h5read(archs4.human, "meta/genes/symbol")

# Extract expression data from ARCHS4
expression <- h5read(archs4.human, "data/expression",
                     index=list(my.sample.locations, 1:length(genes)))
# transpose to get genes as rows and samples as columns
expression <- t(expression)

rownames(expression) <- genes
colnames(expression) <- all.samples.human[my.sample.locations]
colSums(expression) #this shows the sequencing depth for each of the samples you've extracted
archs4.dgelist <- DGEList(expression)
save(archs4.dgelist, file = "archs4.DGEList")
archs4.cpm <- cpm(archs4.dgelist)
colSums(archs4.cpm)

# Filter and normalize the extracted data ----
table(rowSums(archs4.dgelist$counts==0)==147)
keepers <- rowSums(archs4.cpm>1)>=2
archs4.dgelist.filtered <- archs4.dgelist[keepers,]
dim(archs4.dgelist.filtered)
archs4.dgelist.filtered.norm <- calcNormFactors(archs4.dgelist.filtered, method = "TMM")

archs4.filtered.norm.log2.cpm <- cpm(archs4.dgelist.filtered.norm, log=TRUE)

# PART 1: Get sample metadata from ARCHS4 to create a study design file ----
# You don't need to run this code, just showing you what I did to generate the studyDesign_lab11.txt file that provided for this lab
# extract the sample source
sample_source_name <- h5read(archs4.human, "meta/samples/source_name_ch1")
# extract sample title
sample_title <- h5read(archs4.human, name="meta/samples/title")
# extract sample characteristics
sample_characteristics<- h5read(archs4.human, name="meta/samples/characteristics_ch1")

# let's try putting this all together in a study design file
studyDesign <- tibble(Sample_title = sample_title[my.sample.locations],
                      Sample_characteristics = sample_characteristics[my.sample.locations])

studyDesign <- mutate_all(studyDesign, as.vector)

write_tsv(studyDesign,"studyDesign_draft.txt")
# open the file you just created and annotate the samples with the experimental variables (will be obvious from the 'Sample_characteristics' column)

# read this file back in
studyDesign <- read_tsv("studyDesign_lab11.txt")

# PART 2: Principal component analysis (PCA) -------------
#capture experimental variables as factors from this study design
condition <- factor(studyDesign$patient_condition)
skin_type <- factor(studyDesign$skin_type)
condition_skin <- factor(paste(condition,skin_type,sep="_"))
sampleName <- studyDesign$Sample_title

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

pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color=condition_skin) +
  geom_point(size=4) +
  # geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) +
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

# PART 3: GSEA with CAMERA ----
# Set up your experimental design matrix
design <- model.matrix(~0 + condition_skin)
colnames(design) <- levels(condition_skin)

v.archs4.dgelist.filtered.norm <- voom(archs4.dgelist.filtered.norm, design, plot = TRUE)
fit <- lmFit(v.archs4.dgelist.filtered.norm, design)

contrast.matrix <- makeContrasts(psoriasis = psoriasis_lesion - healthy_control,
                                 atopic_dermatitis = atopic_dermatitis_lesion - healthy_control,
                                 psoriasis_vs_AD = psoriasis_lesion - atopic_dermatitis_lesion,
                                 levels=design)


# now test for enrichment using CAMERA
# first, we need to load the gene sets we want to use for GSVA
C2CP <- getGmt("/Users/danielbeiting/Dropbox/MSigDB/c2.cp.v2023.2.Hs.symbols.gmt", geneIdType=SymbolIdentifier())
#extract as a list
C2CP <- geneIds(C2CP)
v.archs4.dgelist.filtered.norm <- voom(archs4.dgelist.filtered.norm, design, plot = TRUE)

camera.res <- camera(v.archs4.dgelist.filtered.norm$E, C2CP, design, contrast.matrix[,2])
camera.df <- as_tibble(camera.res, rownames = "setName")
camera.df


# filter based on FDR and display as interactive table
camera.df <- dplyr::filter(camera.df, FDR<=0.05)

datatable(camera.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2,4,5), digits=2)

#as before, add a variable that maps up/down regulated pathways with phenotype
camera.df <- camera.df %>%
  mutate(phenotype = case_when(
    Direction == "Up" ~ "disease",
    Direction == "Down" ~ "healthy"))

# graph camera results as bubble chart
ggplot(camera.df[1:50,], aes(x=phenotype, y=setName)) +
  geom_point(aes(size=NGenes, color = Direction, alpha=-log10(FDR))) +
  theme_bw()

# PART 4: (BONUS) ----
# what pathways define chronic atopic dermatitis?
# just change your contrast for camera above to 3 and you're good to go!


