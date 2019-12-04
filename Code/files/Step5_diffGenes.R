# Introduction to this script -----------
#the goal of this script is to identify differentially expressed genes (DEGs) or transcripts (DETs)
#you should already know which pairwise comparisons are most important to you
#whether you look for differential expression at the gene or transcript level depends on how you read the Kallisto output into R using TxImport back in Step 1
#if you have no biological replicates, you will NOT be able to leverage statistical tools for differential expression analysis
#Instead, you will ONLY rely on fold changes, and can use the dplyr 'verbs' we discussed in Step 3 to identify genes based on fold change only

# Load packages -----
library(tidyverse) # you know it well by now!
library(limma) # powerful package for differential gene expression using linear modeling
library(edgeR) # well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(gt) # create publication-quality tables from your dataframes
library(DT) # for making interactive tables
library(plotly) # for making interactive plots

# Set up your design matrix ----
treatment <- factor(targets$treatment)
design <- model.matrix(~0 + treatment)
colnames(design) <- levels(treatment)

# NOTE: if you need a paired analysis (a.k.a.'blocking' design) or have a batch effect, the following design is useful
# design <- model.matrix(~block + treatment) 
# this is just an example. 'block' and 'treatment' would need to be objects in your environment

# Model mean-variance trend and fit linear model to data ----
# Use VOOM function from Limma package to model the mean-variance relationship
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)
# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design)

# Contrast matrix ----
#how do cells respond to infection with Crypto?
contrast.matrix <- makeContrasts(infection_with_WT = crypto.wt - uninfected,
                                 infection_with_Mut = crypto.mut - uninfected,
                                 levels=design)

# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)
#stats <- write.fit(ebFit)

# TopTable to view DEGs -----
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=50, sort.by="logFC")

# convert to a tibble
myTopHits <- as_tibble(myTopHits, rownames = "geneSymbol")  

gt(myTopHits)
# TopTable (from Limma) outputs a few different stats:
# logFC, AveExpr, and P.Value should be self-explanatory
# adj.P.Val is your adjusted P value, also known as an FDR (if BH method was used for multiple testing correction)
# B statistic is the log-odds that that gene is differentially expressed. If B = 1.5, then log odds is e^1.5, where e is euler's constant (approx. 2.718).  So, the odds of differential expression os about 4.8 to 1 
# t statistic is ratio of the logFC to the standard error (where the error has been moderated across all genes...because of Bayesian approach)

# Volcano Plots ----
# in topTable function above, set 'number=40000' to capture all genes

# now plot
p <- ggplot(myTopHits, aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneSymbol))) +
  geom_point(size=2) +
  ylim(-0.5,12) +
  #geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  #geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  #geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  annotate("rect", xmin = 1, xmax = 10, ymin = -log10(0.01), ymax = 12, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -1, xmax = -10, ymin = -log10(0.01), ymax = 12, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "C. parvum infected vs. naive (HCT-8 cells)",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# how would you make the volcano plot above interactive?
ggplotly(p)

# decideTests to pull out the DEGs and make Venn Diagram ----
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=2)

# take a look at what the results of decideTests looks like
head(results)
summary(results)
vennDiagram(results, include="both")

# retrieve expression data for your DEGs ----
head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0 | results[,2] !=0,]
head(diffGenes)
dim(diffGenes)
#convert your DEGs to a dataframe using as_tibble
diffGenes.df <- as_tibble(diffGenes, rownames = "geneSymbol")

# create interactive tables to display your DEGs ----
datatable(diffGenes.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: DEGs for infected (wt Crypto) vs control',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(1:6), digits=3)

#write your DEGs to a file
write_tsv(diffGenes.df,"DiffGenes.txt") #NOTE: this .txt file can be directly used for input into Clust (https://github.com/BaselAbujamous/clust)

# the essentials ----
library(tidyverse)
library(limma)
library(edgeR) 
library(gt)
library(DT) 
library(plotly) 

treatment <- factor(targets$treatment)
design <- model.matrix(~0 + treatment)
colnames(design) <- levels(treatment)

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts(infection_with_WT = crypto.wt - uninfected,
                                 infection_with_Mut = crypto.mut - uninfected,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits <- as_tibble(myTopHits, rownames = "geneSymbol")
vplot <- ggplot(myTopHits, aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneSymbol))) +
  geom_point(size=2) +
  ylim(-0.5,12) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  labs(title="Volcano plot",
       subtitle = "C. parvum infected vs. naive (HCT-8 cells)",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot)

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0 | results[,2] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneSymbol")
datatable(diffGenes.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 2: DEGs for infected (wt Crypto) vs control',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(1:10), digits=2)


