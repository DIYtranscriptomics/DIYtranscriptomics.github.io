# Introduction to this script ----
# now that you've read your transcript-level or gene-level data into R, you're ready to begin working with your data
# recall that your abundance data are TPM, while the counts are read counts mapping to each gene or transcript
# Our goals in this script are to filter and normalize data
# This script also introduces us to ggplot2 for plotting, which we'll use to visually depict the impact of filtering and normalization on our data.

# Load packages -----
library(tidyverse) # already know about this from Step 1 script
library(RColorBrewer) # provides access to color palettes for graphics
library(reshape2) # for reshaping dataframes so they play nice with plotting functions
library(genefilter) # as the package name suggests, it's for filtering genes
library(edgeR) # well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(matrixStats) # let's us easily calculate stats on rows or columns of a data matrix
library(cowplot) # allows you to combine multiple plots in one figure

# Examine your data up to this point ----
myTPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts
colSums(myTPM)
colSums(myCounts)

# capture sample labels from the study design file that you worked with and saved as 'targets' in step 1
targets
sampleLabels <- targets$sample

# Generate summary stats for your data ----
# 1st, calculate summary stats for each transcript or gene, and add these to your data matrix
# then use the base R function 'transform' to modify the data matrix (equivalent of Excel's '=')
myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM)
                         )

#look at what you created
head(myTPM.stats)
#produce a scatter plot of the transformed data
ggplot(myTPM.stats, aes(x=SD, y=MED)) +
  geom_hex(shape=16, size=2)
# Experiment with point shape and size
# Experiment with other plot types (e.g. 'geom_hex' instead of 'geom_point')
# How would these graphs change if you log2 converted the data?

# Make a DGElist from your counts, and plot ----
myDGEList <- DGEList(myCounts)
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
# this is our first time using colors in R
# take a look at the color palettes available to you through RColorBrewer
display.brewer.all()
# you can also view only palettes that are colorblind friendly
display.brewer.all(colorblindFriendly = TRUE)
# now select colors from a single palette
myColors <- brewer.pal(nsamples, "Paired")

# 'coerce' your data matrix to a dataframe so that you can use tidyverse tools on it
log2.cpm.df <- as_tibble(log2.cpm)
log2.cpm.df
# add your sample names to this dataframe (we lost these when we read our data in with tximport)
colnames(log2.cpm.df) <- sampleLabels
# use the reshape2 package to 'melt' your dataframe (from wide to tall)
log2.cpm.df.melt <- melt(log2.cpm.df)
log2.cpm.df.melt <- as_tibble(log2.cpm.df.melt)

ggplot(log2.cpm.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.y = "median", 
               geom = "point", 
               shape = 124, 
               size = 6, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) + #using the Sys.time function from base R to print date/time on graph
  coord_flip()
  # what do you think of the distribution of this data?

# Filter your data ----
#first, take a look at how many genes or transcripts have no read counts at all
table(rowSums(myDGEList$counts==0)==9)
# breaking down the line above is a little tricky.  Let's try:
# 1st - 'myDGEList$counts==0' returns a new 'logical matrix' where each observation (gene) is evaluated (TRUE/FALSE) for each variable (sample) as to whether it has zero counts
# 2nd - passing this logical matrix to 'rowsums' allows you to sum the total number of times an observation was 'TRUE' across all samples
# 3rd - adding the '==9' is a simple way of asking how many of the rowsums equaled 9. In other words, how many genes had 0 counts (TRUE) for all samples
# 4th - passing all this to the 'table' function just provides a handy way to summarize the large logical produced in the previous step

# now set some cut-off to get rid of genes/transcripts with low counts
# again using rowSums to tally up the 'TRUE' results of a simple evaluation
# how many genes had more than 1CPM (TRUE) in at least 3 samples
keepers <- rowSums(cpm>1)>=3
# now use base R's simple subsetting method to filter your DGEList based on the logical produced above
myDGEList.filtered <- myDGEList[keepers,]
dim(myDGEList.filtered)

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered) 
colnames(log2.cpm.filtered.df) <- sampleLabels
log2.cpm.filtered.df.melt <- melt(log2.cpm.filtered.df)
log2.cpm.filtered.df.melt <- as_tibble(log2.cpm.filtered.df.melt) 


ggplot(log2.cpm.filtered.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.y = "median", 
               geom = "point", 
               shape = 124, 
               size = 6, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  coord_flip()

# Normalize your data ----
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
# take a look at this new DGEList object...how has it changed?

# use the 'cpm' function from EdgeR to get counts per million from your normalized data
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm)
colnames(log2.cpm.filtered.norm.df) <- sampleLabels
log2.cpm.filtered.norm.df.melt <- melt(log2.cpm.filtered.norm.df)
log2.cpm.filtered.norm.df.melt <- as_tibble(log2.cpm.filtered.norm.df.melt)

ggplot(log2.cpm.filtered.norm.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.y = "median", 
               geom = "point", 
               shape = 124, 
               size = 6, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  coord_flip()

# what if we wanted to put all three violin plots together?
# go back and assign each plot to a variable (rather than printing to the plots viewer)
# we'll use the 'plot_grid' function from the cowplot package to put these together in a figure
plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)

# the essentials ----
library(RColorBrewer) 
library(reshape2) 
library(genefilter)
library(edgeR) 
library(matrixStats)
library(cowplot)

sampleLabels <- targets$sample
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)
nsamples <- ncol(log2.cpm)
myColors <- brewer.pal(nsamples, "Paired")
log2.cpm.df <- as_tibble(log2.cpm)
colnames(log2.cpm.df) <- sampleLabels
log2.cpm.df.melt <- melt(log2.cpm.df)

p1 <- ggplot(log2.cpm.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.y = "median", 
               geom = "point", 
               shape = 124, 
               size = 6, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  coord_flip()

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=3 #user defined
myDGEList.filtered <- myDGEList[keepers,]
log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered) 
colnames(log2.cpm.filtered.df) <- sampleLabels
log2.cpm.filtered.df.melt <- melt(log2.cpm.filtered.df)
log2.cpm.filtered.df.melt <- as_tibble(log2.cpm.filtered.df.melt) 

p2 <- ggplot(log2.cpm.filtered.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.y = "median", 
               geom = "point", 
               shape = 124, 
               size = 6, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) + 
  coord_flip()

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm)
colnames(log2.cpm.filtered.norm.df) <- sampleLabels
log2.cpm.filtered.norm.df.melt <- melt(log2.cpm.filtered.norm.df)

p4 <- ggplot(log2.cpm.filtered.norm.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.y = "median", 
               geom = "point", 
               shape = 124, 
               size = 6, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  coord_flip()

plot_grid(p1, p2, p3, nrow=1, labels = c('A', 'B', 'C'), label_size = 12)

