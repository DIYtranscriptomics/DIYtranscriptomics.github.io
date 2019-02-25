# Introduction to this script ----
# now that you've read your transcript-level or gene-level data into R, you're ready to begin exploring the structure of the data
# goal of this script is to using multivariate statisical approaches to explore the structure of your data
# recall that your abundance data are TPM, while the counts are read counts mapping to each gene or transcript

# Load packages -----
library(tidyverse) 
library(paletteer) # provides access to color palettes for graphics, see https://github.com/EmilHvitfeldt/r-color-palettes for the full list of color palettes
library(reshape2) # for reshaping dataframes
library(trelliscopejs) # must be downloaded from github
library(genefilter) #as the package name suggests, it's for filtering genes
library(edgeR) # also for differential expression, but we only use for the DGEList object
library(matrixStats) # let's us easily calculate stats on any matrix rows or columns

# Identify variables of interest in study design file ----
targets
groups1 <- targets$treatment
groups2 <- targets$treatment2
groups1 <- factor(groups1)
groups2 <- factor(groups2)
sampleLabels <- targets$sample

# Examine your data up to this point ----
myCPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts
# let's quickly graph both matrices to see what we're dealing with

colSums(myCPM)
colSums(myCounts)

# Take a look at the heteroskedasticity of the data ----
# first, calculate row means and standard deviations for each transcript or gene 
# and add these to your data matrix
myCPM.stats <- transform(myCPM, 
                         SD=rowSds(myCPM), 
                         AVG=rowMeans(myCPM),
                         MED=rowMedians(myCPM)
                         )

myCounts.stats <- transform(myCounts, 
                            SD=rowSds(myCounts), 
                            AVG=rowMeans(myCounts), 
                            MED=rowMedians(myCounts)
                            )

head(myCPM.stats)
#produce a scatter plot of the transformed data
ggplot(myCPM.stats, aes(x=SD, y=MED)) +
  geom_point(shape=1, size=4)
# Experiment with point shape and size
# experiment with geom_hex
# how would these graphs change if you log2 converted the data?

# Make a DGElist from your counts ----
myDGEList <- DGEList(Txi_gene$counts)
# take a look at the DGEList object 
myDGEList
#DEGList objects are a good R data file to consider saving to you working directory
save(myDGEList, file = "myDGEList")
#Saved DGEList objects can be easily shared and loaded into an R environment
load(file = "myDGEList")

# use the 'cpm' function from EdgeR to get counts per million
cpm <- cpm(myDGEList) 
colSums(cpm)
log2.cpm <- cpm(myDGEList, log=TRUE)

# Take a look at the distribution of the Log2 CPM
nsamples <- ncol(log2.cpm)
myColors <- brewer.pal(nsamples, "Paired")

# 'coerce' your data matrix to a dataframe so that you can use tidyverse tools on it
log2.cpm.df <- as_tibble(log2.cpm)
log2.cpm.df
# add your sample names to this dataframe (we lost these when we read our data in with tximport)
colnames(log2.cpm.df) <- sampleLabels
# use the reshape2 package to 'melt' your dataframe (from wide to tall)
log2.cpm.df.melt <- melt(log2.cpm.df)
Log2.cpm.df.melt <- as_tibble(log2.cpm.df.melt)
Log2.cpm.df.melt

ggplot(Log2.cpm.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, show.legend = TRUE) +
  stat_summary(fun.y = "median", geom = "point", shape = 95, size = 10, color = "black") +
  theme_bw()
# what do you think of the distribution of this data?

# Filter your data ----
#first, take a look at how many genes or transcripts have no read counts at all
table(rowSums(DGEList$counts==0)==9)

# now set some cut-off to get rid of genes/transcripts with low counts
keepers <- rowSums(cpm>1)>=3
DGEList.filtered <- DGEList[keepers,]
dim(DGEList.filtered)

log2.cpm.filtered <- cpm(DGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered) 
colnames(log2.cpm.filtered.df) <- sampleLabels
log2.cpm.filtered.df.melt <- melt(log2.cpm.filtered.df)
log2.cpm.filtered.df.melt <- as_tibble(log2.cpm.filtered.df.melt)

ggplot(log2.cpm.filtered.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = TRUE, show.legend = TRUE) +
  stat_summary(fun.y = "median", geom = "point", shape = 95, size = 10, color = "black") +
  theme_bw()

# Normalize your data ----
DGEList.filtered.norm <- calcNormFactors(DGEList.filtered, method = "TMM")
# take a look at this new DGEList object...how has it changed?

# use the 'cpm' function from EdgeR to get counts per million from your normalized data
log2.cpm.filtered.norm <- cpm(DGEList.filtered.norm, log=TRUE)

log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm)

colnames(log2.cpm.filtered.norm.df) <- sampleLabels
log2.cpm.filtered.norm.df.melt <- melt(log2.cpm.filtered.norm.df)
log2.cpm.filtered.norm.df.melt <- as_tibble(log2.cpm.filtered.norm.df.melt)

ggplot(log2.cpm.filtered.norm.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = TRUE, show.legend = TRUE) +
  stat_summary(fun.y = "median", geom = "point", shape = 95, size = 10, color = "black") +
  theme_bw()


# the essentials ----
library(paletteer) 
library(reshape2) 
library(genefilter)
library(edgeR) 
library(matrixStats)
library(gridExtra)
groups2 <- targets$treatment2
groups2 <- factor(groups2)
sampleLabels <- targets$sample
myDGEList <- DGEList(Txi_gene$counts)
save(myDGEList, file = "myDGEList")
load(file = "myDGEList")
log2.cpm <- cpm(myDGEList, log=TRUE)
nsamples <- ncol(log2.cpm)
myColors <- brewer.pal(nsamples, "Paired")
log2.cpm.df <- as_tibble(log2.cpm)
colnames(log2.cpm.df) <- sampleLabels
log2.cpm.df.melt <- melt(log2.cpm.df)
Log2.cpm.df.melt <- as_tibble(log2.cpm.df.melt)

p1 <- ggplot(Log2.cpm.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.y = "median", geom = "point", shape = 95, size = 10, color = "black", show.legend = FALSE) +
  labs(y="log2 expression", x = "sample") +
  coord_flip() +
  theme_bw()

keepers <- rowSums(cpm>1)>=3 #user defined
DGEList.filtered <- DGEList[keepers,]
DGEList.filtered.norm <- calcNormFactors(DGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(DGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm)
colnames(log2.cpm.filtered.norm.df) <- sampleLabels
log2.cpm.filtered.norm.df.melt <- melt(log2.cpm.filtered.norm.df)
log2.cpm.filtered.norm.df.melt <- as_tibble(log2.cpm.filtered.norm.df.melt)

p2 <- ggplot(log2.cpm.filtered.norm.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = TRUE, show.legend = FALSE) +
  stat_summary(fun.y = "median", geom = "point", shape = 95, size = 10, color = "black", show.legend = FALSE) +
  labs(y="log2 expression", x = "sample") +
  coord_flip() +
  theme_bw()

grid.arrange(p1, p2, nrow = 1)



