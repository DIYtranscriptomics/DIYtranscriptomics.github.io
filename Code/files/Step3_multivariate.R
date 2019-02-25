# Introduction to this script -----------
# this script walks thorough some basic data wrangling for dealing with any dataframe
# we'll also continue to create publication quality graphics
# we'll start the script with your abundance data from the 

# Load packages ------
library(tidyverse)
library(ggplot2) #graphing package that employs a 'grammar of graphics' approach
library(reshape2) #data manipulation package
library(dplyr) #data manipulation package
library(DT)
library(gt) # for making publication quality tables
library(plotly)

# Data -------
# for this part of the class you'll use your normalized and filtered data in log2 cpm
# make sure you have this object already in your work environment
# if you don't, go back to the Step2 script and generate it
head(log2.cpm.filtered.norm)
#give your data some more informative column names
colnames(log2.cpm.filtered.norm) <- sampleLabels
# Now we need to convert our datamatrix to a dataframe, while preserving the rownames as a new column in this dataframe
mydata.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneSymbol")
write_tsv(mydata.df, "normData.txt") # Note: this is the data you would use as input .gct file for GSEA analysis

# Hierarchical clustering ---------------
#hierarchical clustering can only work on a data matrix, not a data frame
#try using filtered and unfiltered data...how does this change the results
distance <- dist(t(log2.cpm.filtered.norm), method="manhattan") #other dist methods are "maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "complete") #other methods are ward.D, ward.D2, single, complete, average
plot(clusters, labels=sampleLabels)

# Pricipal component analysis (PCA) -------------
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
ggplot(pca.res.df, aes(x=PC1, y=PC2, color=groups1)) +
  geom_point(size=5) +
  theme(legend.position="right") 
#what other variables could you 'paint' onto this PCA plot
#how would this PCA look if you used raw counts (myCounts) instead of log2 CPM?

# Create a PCA 'small multiples' chart ----
# this is another way to view PCA laodings to understand impact of each sample on each pricipal component
melted <- cbind(groups1, melt(pca.res$x[,1:4]))
head(melted)
#look at your 'melted' data
ggplot(melted) +
  geom_bar(aes(x=Var1, y=value, fill=groups1), stat="identity") +
  facet_wrap(~Var2) +
  theme(axis.text.x = element_text(angle = 90))
#facet_trelliscope(~Var2) #add this trelliscope layer to the ggplot to have some fun


# use dplyr 'mutate' function to add new columns based on existing data -------
mydata.df <- mutate(mydata.df,
                   uninfected.AVG = (uninf_rep1 + uninf_rep2 + uninf_rep1)/3, 
                   crypto.wt.AVG = (crypto.wt_rep1 + crypto.wt_rep2 + crypto.wt_rep3)/3,
                   crypto.mut.AVG = (crypto.mut_rep1 + crypto.mut_rep2 + crypto.mut_rep3)/3,
                   #now make columns comparing each of the averages above that you're interested in
                   LogFC.crypto.wt_vs_uninfected = (crypto.wt.AVG - uninfected.AVG),
                   LogFC.crypto.mut_vs_uninfected = (crypto.mut.AVG - uninfected.AVG)) %>%
  mutate_if(is.numeric, ~round(., 2))
#why is this type of approach to managing a spreadsheet useful?

#now look at this modified data table
mydata.df

# using dplyr 'arrange' and 'select' to sort your dataframe based on any variable -----
# first, we'll use dplyr "arrange" function to sort rows based on the values in a column of interest
# then we'll display 'select' only the columns we're interested in seeing
mydata.sort <- mydata.df %>%
  dplyr::arrange(desc(LogFC.crypto.wt_vs_uninfected)) %>% #note that this is the first time you've seen the 'pipe' operator
  dplyr::select(geneSymbol, LogFC.crypto.wt_vs_uninfected)


# use dplyr "filter" and "select" functions to pick out genes of interest  ----
#ways to tweek the 'select' function
#use : between two column names to select all columns between
#use 'contains', 'starts_with' or 'ends_with' to modify how you select
#can refer to columns using exact name or numerical indicator
#use boolean operators such as '&' (and), '|' (or), '==' (equal to), '!' (not)
mydata.filter <- mydata.df %>%
  dplyr::filter(geneSymbol=="MX1" | geneSymbol=="IRF1" | geneSymbol=="OAS2" | geneSymbol=="IRF3") %>%
  dplyr::select(geneSymbol, uninfected.AVG, crypto.wt.AVG)
mydata.filter

# you can also filter based on any regular expression
mydata.grep <- mydata.df %>%
  dplyr::filter(grepl('CXCL|IFI', geneSymbol)) %>%
  dplyr::select(geneSymbol, uninfected.AVG, crypto.wt.AVG, LogFC.crypto.wt_vs_uninfected) %>%
  dplyr::arrange(desc(geneSymbol))

mydata.grep
#let's convert this tibble to a publication quality table using the gt package
gt(mydata.grep)
#now with a few more options
gt(mydata.grep %>%
     fmt_number(columns = vars(crypto.wt.AVG), decimals = 1) %>%
     tab_header(
       title = md("**My favorite genes!**"),
       subtitle = md("but just *some* of them")
     ))

# make an interactive table ----
datatable(mydata.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'my cool table)',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(1:6), digits=3)

# make a simple scatter plot -----
ggplot(mydata.df, aes(x=uninfected.AVG, y=crypto.wt.AVG)) +
  geom_point(shape=3) +
  geom_point(size=1)

# the default ggplot theme is not very aesthetically pleasing and there are lots of ways to customize the way your plot looks
# perhaps one of the easiest ways to simplify/clean the aesthetics for publicaiton is to use the cowplot package
# cowplot, once loaded, will change the default theme for all ggplots your produce
library(cowplot) 
# see what happens as we regenerate the previous plot
ggplot(mydata.df, aes(x=uninfected.AVG, y=crypto.wt.AVG)) +
  geom_point(shape=3) +
  geom_point(size=1)


# now we'll use the plotly package to convert our static scatter plot to an interactive
# begin by storing your ggplot object in a variable
myplot <- ggplot(mydata.df, aes(x=uninfected.AVG, y=crypto.wt.AVG)) +
  geom_point(shape=3) +
  geom_point(size=1)
#now use the ggplotly function from the plotly package to convert this ggplot object into an interactive
ggplotly(myplot)

#let's customize this graphic by adding a more informative tooltip
myplot <- ggplot(mydata.df, aes(x=uninfected.AVG, y=crypto.wt.AVG, text = paste("Symbol:", geneSymbol))) +
  geom_point(shape=3) +
  geom_point(size=1)

ggplotly(myplot)

# the essentials ----
library(tidyverse)
library(ggplot2) 
library(reshape2) 
library(dplyr) 
library(DT)
library(gt) 
library(plotly)

mydata.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneSymbol")
write_tsv(mydata.df, "normData.txt")
#user defined
mydata.df <- mutate(mydata.df,
                    uninfected.AVG = (uninf_rep1 + uninf_rep2 + uninf_rep1)/3, 
                    crypto.wt.AVG = (crypto.wt_rep1 + crypto.wt_rep2 + crypto.wt_rep3)/3,
                    crypto.mut.AVG = (crypto.mut_rep1 + crypto.mut_rep2 + crypto.mut_rep3)/3,
                    #now make columns comparing each of the averages above that you're interested in
                    LogFC.crypto.wt_vs_uninfected = (crypto.wt.AVG - uninfected.AVG),
                    LogFC.crypto.mut_vs_uninfected = (crypto.mut.AVG - uninfected.AVG)) %>%
  mutate_if(is.numeric, ~round(., 2))

datatable(mydata.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'my cool table)',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(1:6), digits=3)

pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
x <- pca.res$rotation 
pc.var<-pca.res$sdev^2
pc.per<-round(pc.var/sum(pc.var)*100, 1)

pca.res.df <- as_tibble(pca.res$x)
p3 <- ggplot(pca.res.df, aes(x=PC1, y=PC2, color=groups1)) +
  geom_point(size=5) +
  theme(legend.position="right") 
ggplotly(p3)

