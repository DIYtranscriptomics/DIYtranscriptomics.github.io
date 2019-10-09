# Introduction to this script -----------
# this script walks thorough some basic data wrangling for dealing with any dataframe
# we'll also continue to create publication -uality graphics
# This script starts with your filtered and normalized abundance data from the Step 2 script.

# Load packages ------
library(tidyverse)
library(reshape2)
library(DT) # for making interactive tables
library(plotly) # for making interactive plots
library(skimr) # great tool for summarizing dataframes
library(gt) # create publication-quality tables from your dataframes

# Identify variables of interest in study design file ----
targets
groups1 <- targets$treatment
groups2 <- targets$treatment2
treatment <- factor(groups2)
sampleLabels <- targets$sample

# Prepare your data -------
# for this part of the class you'll use your normalized and filtered data in log2 cpm
# make sure you have this object already in your work environment
# if you don't, go back to the Step2 script and generate it
head(log2.cpm.filtered.norm)
# Now we need to convert our datamatrix to a dataframe, while preserving the rownames as a new column in this dataframe
mydata.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneSymbol")
colnames(mydata.df) <- c("geneSymbol", sampleLabels)
skim(mydata.df)

# Hierarchical clustering ---------------
#hierarchical clustering can only work on a data matrix, not a data frame
#try using filtered and unfiltered data...how does this change the result?
#try other distance methods (e.g. switch from 'maximum' to 'euclidean')...how does this change the result?
distance <- dist(t(log2.cpm.filtered.norm), method="maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(clusters, labels=sampleLabels)

# Principal component analysis (PCA) -------------
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
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
ggplot(pca.res.df, aes(x=PC1, y=PC2, color=targets$treatment)) +
  geom_point(size=4) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

#what other variables could you 'paint' onto this PCA plot
#how would this PCA look if you used raw counts (myCounts) instead of log2 CPM?
#what are the disadvantages of looking at a PCA result using such a simple XY plot?

# Create a PCA 'small multiples' chart ----
# this is another way to view PCA laodings to understand impact of each sample on each pricipal component
melted <- cbind(groups1, melt(pca.res$x[,1:4]))
head(melted) #look at your 'melted' data
colnames(melted) <- c('group', 'treatment', 'PC', 'loadings')
ggplot(melted) +
  geom_bar(aes(x=treatment, y=loadings, fill=group), stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_flip() +
  theme_bw()

# Use dplyr 'mutate' function to add new columns based on existing data -------
mydata.df <- mutate(mydata.df,
                   uninfected.AVG = (uninf_rep1 + uninf_rep2 + uninf_rep1)/3, 
                   crypto.wt.AVG = (crypto.wt_rep1 + crypto.wt_rep2 + crypto.wt_rep3)/3,
                   crypto.mut.AVG = (crypto.mut_rep1 + crypto.mut_rep2 + crypto.mut_rep3)/3,
                   #now make columns comparing each of the averages above that you're interested in
                   LogFC.crypto.wt_vs_uninfected = (crypto.wt.AVG - uninfected.AVG),
                   LogFC.crypto.mut_vs_uninfected = (crypto.mut.AVG - uninfected.AVG)) %>%
  mutate_if(is.numeric, round, 2)

#now look at this modified data table
mydata.df

# Use dplyr 'arrange' and 'select' to sort your dataframe based on any variable -----
# first, we'll use dplyr "arrange" function to sort rows based on the values in a column of interest
# then we'll display 'select' only the columns we're interested in seeing
mydata.sort <- mydata.df %>%
  dplyr::arrange(desc(LogFC.crypto.wt_vs_uninfected)) %>% #note that this is the first time you've seen the 'pipe' operator
  dplyr::select(geneSymbol, LogFC.crypto.wt_vs_uninfected)

# Use dplyr "filter" and "select" functions to pick out genes of interest  ----
#ways to tweek the 'select' function
#use ':' between two column names to select all columns between
#use 'contains', 'starts_with' or 'ends_with' to modify how you select
#can refer to columns using exact name or numerical indicator
#use boolean operators such as '&' (and), '|' (or), '==' (equal to), '!' (not)
mydata.filter <- mydata.df %>%
  dplyr::filter(geneSymbol=="MX1" | geneSymbol=="IRF1" | geneSymbol=="OAS2" | geneSymbol=="IRF3") %>%
  dplyr::select(geneSymbol, uninfected.AVG, crypto.wt.AVG, LogFC.crypto.wt_vs_uninfected)

# you can also filter based on any regular expression
mydata.grep <- mydata.df %>%
  dplyr::filter(grepl('CXCL|IFI', geneSymbol)) %>%
  dplyr::select(geneSymbol, uninfected.AVG, crypto.wt.AVG, LogFC.crypto.wt_vs_uninfected) %>%
  dplyr::arrange(desc(geneSymbol))

#let's convert this tibble to a publication quality table using the gt package
gt(mydata.grep)
#now with a few more options
gt(mydata.grep %>%
     fmt_number(columns = vars(crypto.wt.AVG), decimals = 1) %>%
     tab_header(
       title = md("**My favorite genes!**"),
       subtitle = md("but just *some* of them")))

# Make an interactive table ----
datatable(mydata.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'my cool table)',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100")))

# Make an interactive scatter plot -----
ggplot(mydata.df, aes(x=uninfected.AVG, y=crypto.wt.AVG)) +
  geom_point(shape=16) +
  geom_point(size=1) +
  theme_bw()

# now we'll use the plotly package to convert our static scatter plot to an interactive
# begin by storing your ggplot object
myplot <- ggplot(mydata.df, aes(x=uninfected.AVG, y=crypto.wt.AVG)) +
  geom_point(shape=16) +
  geom_point(size=1) +
  ggtitle("Infected vs. naive") +
  theme_bw()

#now use the ggplotly function from the plotly package to convert this ggplot object into an interactive
ggplotly(myplot)

#let's customize this graphic by adding a more informative tooltip
myplot <- ggplot(mydata.df, aes(x=uninfected.AVG, y=crypto.wt.AVG,
                                text = paste("Symbol:", geneSymbol))) +
  geom_point(shape=16) +
  geom_point(size=1) +
  ggtitle("Infected vs. naive") +
  theme_bw()

ggplotly(myplot)

# the essentials ----
library(tidyverse)
library(reshape2)
library(DT)
library(gt)
library(plotly)

mydata.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneSymbol")
colnames(mydata.df) <- c("geneSymbol", sampleLabels)

mydata.df <- mutate(mydata.df,
                    uninfected.AVG = (uninf_rep1 + uninf_rep2 + uninf_rep1)/3, 
                    crypto.wt.AVG = (crypto.wt_rep1 + crypto.wt_rep2 + crypto.wt_rep3)/3,
                    crypto.mut.AVG = (crypto.mut_rep1 + crypto.mut_rep2 + crypto.mut_rep3)/3,
                    #now make columns comparing each of the averages above that you're interested in
                    LogFC.crypto.wt_vs_uninfected = (crypto.wt.AVG - uninfected.AVG),
                    LogFC.crypto.mut_vs_uninfected = (crypto.mut.AVG - uninfected.AVG)) %>%
  mutate_if(is.numeric, round, 2)

datatable(mydata.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'my cool table)',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100")))

distance <- dist(t(log2.cpm.filtered.norm), method="maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(clusters, labels=sampleLabels)

pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
x <- pca.res$rotation 
pc.var<-pca.res$sdev^2
pc.per<-round(pc.var/sum(pc.var)*100, 1)

pca.res.df <- as_tibble(pca.res$x)

pca.plot <- ggplot(pca.res.df, aes(x=PC1, y=PC2, color=groups1)) +
  geom_point(size=5) +
  theme(legend.position="right") 
ggplotly(pca.plot)
