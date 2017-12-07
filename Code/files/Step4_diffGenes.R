# Introduction to this script -----------
#the goal of this script is to identify differentially expressed genes (DEGs) or transcripts (DETs)
#you should already know what pairwise comparisons are most important to you
#whether you look for differential expression at the gene or transcript level depends on how you read the Kallisto output into R using TxImport
#if you have no biological replicates, you will NOT be able to leverage statistical tools for differential expression analysis
#Instead, you will ONLY rely on fold changes, and can use the dplyr 'verbs' we discussed in the last class to identify genes based soley on fold changes

# Load packages -----
library(limma) #powerful package for differential gene expression using linear modeling
library(edgeR)
library(sleuth)


# OPTION 1: DE analysis using Sleuth ----
# Make sure to do the following first:
# 1. make sure your studyDesign file has the correct path to each Kallisto output folder
# 2. change your design matrix to set an intercept
# 3. if you used Salmon or Sailfish, instead of Kallisto, for read mapping, you can still use Sleuth, but you'll need the Wasabi R package
targets <- read.table("Crypto_studyDesign.txt", row.names=NULL, header = T, as.is = T)
design <- model.matrix(~groups)

# Now you're ready to construct a sleuth object
mySleuth <- sleuth_prep(targets, 
                        ~ treatment, 
                        target_mapping = Tx,
                        read_bootstrap_tpm=TRUE, 
                        extra_bootstrap_summary=TRUE) 

# fit a linear model to the data
so <- sleuth_fit(mySleuth, design)
# use a wald test (WT) to identify genes that are differentially expressed
so.WT_crypto <- sleuth_wt(so, which_beta = "groupswt_crypto")
sleuth_live(so.WT_crypto)

# OPTION 2: DE analysis using Limma/VOOM (alternatively, EdgeR or DESeq2) -----
# first create a DGEList object from your original count data using the DGEList function from EdgeR
DGEList.filtered.norm


#set up your design matrix
groups 
#groups <- relevel(groups, "control") #may need to use 'relevel' function
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

#normalize your data using the mean-variance relationship using the VOOM function from Limma
v.DEGList.filtered.norm <- voom(DGEList.filtered.norm, design, plot = TRUE)
# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design)

# Contrast matrix ----
#how do cells respond to infection with Crypto?
contrast.matrix <- makeContrasts(WT_crypto_infection = wt_crypto - control,
                                 trans_crypto_infection = trans_crypto - control,
                                 levels=design)

# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)
#stats <- write.fit(ebFit)

# TopTable to view DEGs -----
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=2000, sort.by="logFC")
head(myTopHits)
# Volcano Plots ----------
# Make a basic volcano plot
with(myTopHits, plot(logFC, -log10(adj.P.Val), pch=20, main="Volcano plot"))
with(subset(myTopHits, adj.P.Val<.05 & abs(logFC)>1), points(logFC, -log10(adj.P.Val), pch=20, col="red"))

# make an interactive volcano plot 
#first, move rownames into the dataframe
library(tibble)
myTopHits <- rownames_to_column(myTopHits, "geneID")
head(myTopHits)

#set-up your tool tip
tooltip <- function(data, ...) {
  paste0("<b>","Symbol: ", data$geneID, "</b><br>",
         "LogFC: ", data$logFC, "<br>",
         "FDR: ", data$adj.P.Val)
}


#plot the interactive graphic
library(ggvis)
myTopHits %>% 
  ggvis(x= ~logFC, y= ~-log10(adj.P.Val), key := ~geneID) %>% 
  add_tooltip(tooltip)

#make some beautiful tables to display your DEGs
library(DT)
datatable(myTopHits, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: DEGs for infected (wt Crypto) vs control',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(1:6), digits=3)

# decideTests to pull out the DEGs and make Venn Diagram ----
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=1)

# take a look at what the results of decideTests looks like
tail(results)
summary(results)
vennDiagram(results, include="up")

# retrieve expression data for your DEGs ----
head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0 | results[,2] !=0,]
head(diffGenes)
dim(diffGenes)
#write your DEGs to a file
write.csv(diffGenes,"DiffGenes.csv")
write.table(diffGenes,"DiffGenes.txt", sep = "\t")

# OPTIONAL: paired design and correcting for known batch effects ----
# if you need a paired analysis (a.k.a.'blocking' design)
# design <- model.matrix(~block+treatment) #this is just an example. 'block' and 'treatment' would need to be objects in your environment
# if your exploratory analysis showed a batch effect
library(sva)
batch <- targets$batch
normData.batchCorrected <- ComBat(normData, batch, design)
#now use this batch corrected data moving forward


