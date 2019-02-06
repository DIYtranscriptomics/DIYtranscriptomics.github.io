# Introduction to this script -----------
#the goal of this script is to identify differentially expressed genes (DEGs) or transcripts (DETs)
#you should already know which pairwise comparisons are most important to you
#whether you look for differential expression at the gene or transcript level depends on how you read the Kallisto output into R using TxImport back in Step 1
#if you have no biological replicates, you will NOT be able to leverage statistical tools for differential expression analysis
#Instead, you will ONLY rely on fold changes, and can use the dplyr 'verbs' we discussed in Step 3 to identify genes based on fold change only

# Load packages -----
library(tidyverse)
library(limma) #powerful package for differential gene expression using linear modeling
library(edgeR) #another great package for differential gene expression analysis
library(scatterD3) #creates interactive plots using ggplot commands
library(DT) #creates interactive datatables
library(sva) #includes many functions for handling sources of variance (e.g. batch effects)

# Set up your design matrix ----
# remember how we read in our study design and capture our treatment variable (from Step 1 script)
targets <- read_tsv("Crypto_studyDesign.txt")
targets
groups1 <- targets$treatment
groups2 <- targets$treatment2
groups1 <- factor(groups1)
groups2 <- factor(groups2)
groups1 <- relevel(groups1, "uninfected") #may need to use 'relevel' function
design <- model.matrix(~0 + groups1)
colnames(design) <- levels(groups1)

# Model mean-variance trend and fit linear model to data ----
# Use VOOM function from Limma package to model the mean-variance relationship
v.DEGList.filtered.norm <- voom(DGEList.filtered.norm, design, plot = TRUE)
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
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=10, sort.by="logFC")
myTopHits

# Volcano Plots ----
# first, move rownames into the dataframe
myTopHits <- rownames_to_column(myTopHits, "geneID")
myTopHits

# now plot
ggplot(myTopHits, aes(y=-log10(adj.P.Val), x=logFC)) +
  geom_point(size=3) +
  theme_linedraw() +
  theme(legend.position="right", axis.title = element_text(size = 19), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.text = element_text(size = 15),
        #axis.title.x=element_blank(),
        axis.text.x = element_text(size=14, colour = "black"),
        axis.text.y = element_text(size=14, colour = "black")) +
  ylim(-0.5,15) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1)
#annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 14, alpha=.2, fill="#BE684D") +
#annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 14, alpha=.2, fill="#2C467A")


# make an interactive volcano plot 
myTopHits <- mutate(myTopHits, 
                    log10Pval = -log10(adj.P.Val), 
                    adj.P.Val = round(adj.P.Val, 2),
                    logFC = round(logFC, 2),
                    geneID = geneID)

tooltip1 <- paste("<b>","Symbol: ", myTopHits$geneID, "</b><br>",
                  "<b>","LogFC: ","</b>", myTopHits$logFC, "<br>",
                  "<b>","FDR: ","</b>", myTopHits$adj.P.Val, "<br>")

scatterD3(myTopHits, x = logFC, y = log10Pval,
          lasso = TRUE,
          xlab = "logFC", 
          ylab = "-log10(adj.P.val)",
          colors = "#C9B3D6",
          point_opacity = 0.7,
          caption = list(title = "Volcano plot - infected vs. naive"),
          tooltip_text = tooltip1, hover_size = 3)

# create interactive tables to display your DEGs ----
datatable(myTopHits, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: DEGs for infected (wt Crypto) vs control',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(1:6), digits=3)

# decideTests to pull out the DEGs and make Venn Diagram ----
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=2)

# take a look at what the results of decideTests looks like
tail(results)
summary(results)
vennDiagram(results, include="both")

# retrieve expression data for your DEGs ----
head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0 | results[,2] !=0,]
head(diffGenes)
dim(diffGenes)
#convert your DEGs to a dataframe using as_tibble
diffGenes.df <- as_as_tibble(diffGenes, rownames = "geneSymbol")

#write your DEGs to a file
write_tsv(diffGenes.df,"DiffGenesTEST.txt") #NOTE: this .txt file can be directly used for input into Clust (https://github.com/BaselAbujamous/clust)

# OPTIONAL: other study designs ----

# if you need a paired analysis (a.k.a.'blocking' design)
design <- model.matrix(~block + treatment) #this is just an example. 'block' and 'treatment' would need to be objects in your environment

# if your exploratory analysis showed a batch effect
library(sva)
batch <- targets$batch
normData.batchCorrected <- ComBat(normData, batch, design)


