# Introduction to this script -----------
# this script walks thorough techniques for data exploration and expands on last week's data wrangling theme
# we'll also continue to create publication-quality graphics
# This script starts with your filtered and normalized abundance data from the Step 2 script.

# Load packages ------
library(tidyverse) # you're familiar with this fromt the past two lectures
library(DT) # for making interactive tables
library(plotly) # for making interactive plots
library(gt) # A layered 'grammar of tables' - think ggplot, but for tables

# Identify variables of interest in study design file ----
targets
group <- targets$group
group <- factor(group)

# Prepare your data -------
# for this part of the class you'll use your normalized and filtered data in log2 cpm
# make sure you have this object already in your work environment
# if you don't, go back to the Step2 script and generate it
log2.cpm.filtered.norm.df

# Hierarchical clustering ---------------
#hierarchical clustering can only work on a data matrix, not a data frame
#try using filtered and unfiltered data...how does this change the result?
#try other distance methods (e.g. switch from 'maximum' to 'euclidean')...how does this change the result?
distance <- dist(t(log2.cpm.filtered.norm), method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(clusters, labels=sampleLabels)

# Principal component analysis (PCA) -------------
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
#look at the PCA result (pca.res) that you just created
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
pca.res$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
pca.res$x # 'x' shows you how much each sample influenced each PC (called 'loadings')
#note that these have a magnitude and a direction (this is the basis for making a PCA plot)
screeplot(pca.res) # A screeplot is a standard way to view eigenvalues for each PCA
pc.var<-pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc.per

# Visualize your PCA result ------------------
#lets first plot any two PCs against each other
#We know how much each sample contributes to each PC (loadings), so let's plot
pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels) +
  geom_point(size=4) +
  # geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  # coord_fixed() +
  theme_bw()

# Let's discuss and iteratively refine the PCA code and plot from above
# First, take note of the fact that we can use information from our PCA analysis to label our axes
# Remember that PCA is unsupervised, so knows nothing about group assignment (healthy vs disease)
# But *we* know, and so we can use this knowledge to enhance the plot.  Add a 'color=group' mapping to the aes of the plot above
# Can we figure out the identity of the outlier?  We have already provided samplelabel mapping in aes, so just uncomment the 'geom_label()'
# Uncomment 'coord_fixed()' to apply the correct aspect ratio
# Uncomment 'stat_ellipse()' to see how you can circle clusters on the PCA
# How would this PCA look if you used raw counts (myCounts) instead of log2 CPM?
# What are the disadvantages of looking at a PCA result using such a simple XY plot?

# Create a PCA 'small multiples' chart ----
# this is another way to view PCA laodings to understand impact of each sample on each pricipal component
pca.res.df <- pca.res$x[,1:4] %>% # note that this is the first time you've seen the 'pipe' operator from the magrittr package
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group)
  
pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=group) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()



# Use dplyr 'verbs' to modify our dataframe ----
# use dplyr 'mutate' function to add new columns based on existing data
mydata.df <- log2.cpm.filtered.norm.df %>% 
  mutate(healthy.AVG = (HS01 + HS02 + HS03 + HS04 + HS05)/5,
         disease.AVG = (CL08 + CL10 + CL11 + CL12 + CL13)/5,
         #now make columns comparing each of the averages above that you're interested in
         LogFC = (disease.AVG - healthy.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

#now look at this modified data table
mydata.df

# Use dplyr 'arrange' and 'select' to sort your dataframe based on any variable
# first, we'll use dplyr "arrange" function to sort rows based on the values in a column of interest
# then we'll display 'select' only the columns we're interested in seeing
mydata.sort <- mydata.df %>%
  dplyr::arrange(desc(LogFC)) %>% 
  dplyr::select(geneID, LogFC)

# Use dplyr "filter" and "select" functions to pick out genes of interest 
# ways to tweak the 'select' function:
# use ':' between two column names to select all columns between
# use 'contains', 'starts_with' or 'ends_with' to modify how you select
# can refer to columns using exact name or numerical indicator
# use boolean operators such as '&' (and), '|' (or), '==' (equal to), '!' (not)
mydata.filter <- mydata.df %>%
  dplyr::filter(geneID=="MMP1" | geneID=="GZMB" | geneID=="IL1B" | geneID=="GNLY" | geneID=="IFNG"
                | geneID=="CCL4" | geneID=="PRF1" | geneID=="APOBEC3A" | geneID=="UNC13A" ) %>%
  dplyr::select(geneID, healthy.AVG, disease.AVG, LogFC) %>%
  dplyr::arrange(desc(LogFC))

# you can also filter based on any regular expression
mydata.grep <- mydata.df %>%
  dplyr::filter(grepl('CXCL|IFI', geneID)) %>%
  dplyr::select(geneID, healthy.AVG, disease.AVG, LogFC) %>%
  dplyr::arrange(desc(geneID))

# Produce publication-quality tables using the gt package ----
gt(mydata.filter)
# now with a few more options
mydata.filter %>%
  gt() %>%
  fmt_number(columns=2:4, decimals = 1) %>%
  tab_header(title = md("**Regulators of skin pathogenesis**"),
             subtitle = md("*during cutaneous leishmaniasis*")) %>%
  tab_footnote(
    footnote = "Deletion or blockaid ameliorates disease in mice",
    locations = cells_body(
      columns = geneID,
      rows = c(6, 7))) %>% 
  tab_footnote(
    footnote = "Associated with treatment failure in multiple studies",
    locations = cells_body(
      columns = geneID,
      rows = c(2:9))) %>%
  tab_footnote(
    footnote = "Implicated in parasite control",
    locations = cells_body(
      columns = geneID,
      rows = c(2))) %>%
  tab_source_note(
    source_note = md("Reference: Amorim *et al*., (2019). DOI: 10.1126/scitranslmed.aar3619"))


# Make an interactive table using the DT package ----
datatable(mydata.df[,c(1,12:14)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100")))

# Make an interactive scatter plot with plotly -----
# begin by storing your ggplot object
myplot <- ggplot(mydata.df) + 
  aes(x=healthy.AVG, y=disease.AVG) +
  geom_point(shape=16, size=1) +
  ggtitle("disease vs. healthy") +
  theme_bw()

#now use the ggplotly function from the plotly package to convert this ggplot object into an interactive plot
ggplotly(myplot)

#let's customize this graphic by adding a more informative mouseover tooltip
myplot <- ggplot(mydata.df) +
  aes(x=healthy.AVG, y=disease.AVG, 
      text = paste("Symbol:", geneID)) +
  geom_point(shape=16, size=1) +
  ggtitle("disease vs. healthy") +
  theme_bw()

ggplotly(myplot)

# the essentials ----
library(tidyverse)
library(DT)
library(gt)
library(plotly)

group <- targets$group
group <- factor(group)

pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot)

mydata.df <- mutate(log2.cpm.filtered.norm.df,
                    healthy.AVG = (HS01 + HS02 + HS03 + HS04 + HS05)/5, 
                    disease.AVG = (CL08 + CL10 + CL11 + CL12 + CL13)/5,
                    #now make columns comparing each of the averages above that you're interested in
                    LogFC = (disease.AVG - healthy.AVG)) %>% #note that this is the first time you've seen the 'pipe' operator
  mutate_if(is.numeric, round, 2)

datatable(mydata.df[,c(1,12:14)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         #dom = "Blfrtip", 
                         #buttons = c("copy", "csv", "excel"),
                         lengthMenu = c("10", "25", "50", "100")))
