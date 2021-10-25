# Introduction to this script -----------
#this script creates heatmaps from your differentially expressed genes or transcripts
#and selects modules of co-expressed genes based on pearson correlations

# Load packages -----
library(tidyverse)
library(limma) #we only use limma in this script for the 'avearrays' function
library(RColorBrewer) #need colors to make heatmaps
library(gplots) #the heatmap2 function in this package is a primary tool for making heatmaps
library(gameofthrones) #because...why not.  Install using 'devtools::install_github("aljrico/gameofthrones")'
library(heatmaply) #for making interactive heatmaps using plotly
library(d3heatmap) #for making interactive heatmaps using D3

# Choose your color pallette ----
#Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]).
myheatcolors1 <- bluered(75) #this is from the 'colorpanel' function in gplots (same package that heatmap.2 comes from)
display.brewer.all()
display.brewer.all(colorblindFriendly = TRUE)

# lots of easy ways to assemble your own color palette, including:
# 1). use 'colorRampPalette' function from the grDevices package
myheatcolors2 <- colorRampPalette(colors=c("blue","white","red"))(100)
# 2). use rcolorbrewer to choose any palette by name and n colors from that palette
myheatcolors3 <- brewer.pal(name="RdBu", n=11)
# 3). paste in your own hex codes using the Sip app (or other tools)
myheatcolors3 <- c("#fed976", "#268f9c")
# 4). have some fun with outside color packages (e.g. GameOfThrones)
got_palette <- got(75, option = "Arya")

# Data----
# you can make a heatmap out of any datamatrix
# we'll use our 'diffgenes' datamatrix that was produced at the end of the last class in the Step 5 script
# as a reminder, this was produced as follows
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=2)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
dim(diffGenes)

# Cluster DEGs ----
#begin by clustering the genes (rows) in each set of differentially expressed genes
# we use the 'cor' function and the pearson method for finding all pairwise correlations of genes
# '1-cor' converts this to a 0-2 scale for each of these correlations, which can then be used to calculate a distance matrix using 'as.dist'
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") 
# hierarchical clustering is a type of unsupervised clustering. Related methods include K-means, SOM, etc 
# unsupervised methods are blind to sample/group identity
# in contrast, supervised methods 'train' on a set of labeled data.  
# supervised clustering methods include random forest, and artificial neural networks

#now cluster your samples (columns)
#we may not acutally use this clustering result, but it's good to have just in case
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete") #cluster columns by spearman correlation
#note: we use Spearman, instead of Pearson, for clustering samples because it gives equal weight to highly vs lowly expressed transcripts or genes
#see lecture slides for an example using a mock dataset.

#Cut the resulting tree and create color vector for clusters.  
#Vary the cut height (h =) to give more or fewer clusters, or use force k= number of clusters
#we'll look at these clusters in more detail later
module.assign <- cutree(clustRows, k=2)

#now assign a color to each module (makes it easy to identify and manipulate)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 

# Produce a static heatmap of DEGs ----
#plot the hclust results as a heatmap
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors3), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20)) 

#what do the colors represent in this heatmap?
#what happens when you change scale='none'

# Make interactive heatmap ----
#first, we'll make an interactive heatmap using plotly (https://plot.ly/)
heatmaply(diffGenes,
          colors = myheatcolors2,
          Rowv=as.dendrogram(clustRows),
          RowSideColors=module.color,
          #showticklabels=c(FALSE,FALSE),
          scale='row')

# now let's try using D3 to create an html widget version of our heatmap
d3heatmap(diffGenes,
          colors = myheatcolors2,
          Rowv=as.dendrogram(clustRows),
          row_side_colors = module.color,
          scale='row')

# OPTIONAL: simplify heatmap ----
#notice that the heatmap includes ALL the columns from your dataset
#a useful way to simplify heatmaps, especially when there are many conditions, is to average your biological replicates and display only one column per condition
#rerun the heatmap script above using diffData.AVG as input instead of diffData
colnames(diffGenes) <- targets$group

#now an old function from the limma package to average your replicates 
diffGenes.AVG <- avearrays(diffGenes)

##alternatively, decide exactly which columns you want to show, and modify the heatmap accordingly
#this is how it would look using base R
#diffGenes.subset <- diffGenes[,c(1,4,7)]
##now repeat heatmap using only these selected columns

# View modules of co-regulated genes ----
# view your color assignments for the different clusters
names(module.color) <- names(module.assign) 

module.assign.df <- as_tibble(as.list(module.assign))
module.assign.pivot <- pivot_longer(module.assign.df, # dataframe to be pivoted
                          cols = 1:2307, # column names to be stored as a SINGLE variable
                          names_to = "geneID", # name of that new variable (column)
                          values_to = "module") # name of new variable (column) storing all the values (data)

module.assign.pivot <- module.assign.pivot %>%
  mutate(moduleColor = case_when(
    module == 1 ~ "#FF9900",
    module == 2 ~ "#FF0099"))


ggplot(module.assign.pivot) +
  aes(module) +
  geom_bar(aes(fill=moduleColor)) +
  theme_bw()

#choose a cluster(s) of interest by selecting the corresponding number based on the previous graph
modulePick <- 2 #use 'c()' to grab more than one cluster from the heatmap.  e.g., c(1,2)
#now we pull out the genes from this module using a fancy subsetting operation on a named vector
myModule <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub <- hclust(as.dist(1-cor(t(myModule), method="pearson")), method="complete") 

# Create heatmap for chosen sub-cluster.
heatmap.2(myModule, 
          Rowv=as.dendrogram(hrsub), 
          Colv=NA, 
          labRow = NA,
          col=rev(myheatcolors3), scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20)) 

# Export modules for downstream analysis ----
#prints out genes in the order you see them in the cluster
moduleSymbols <- tibble(geneID = rev(hrsub$labels[hrsub$order]))
moduleData <- diffGenes[moduleSymbols$geneID,]
moduleData.df <- as_tibble(moduleData, rownames = "geneSymbol")
write_tsv(moduleData.df,"module_upRegulated.tsv")

# OPTIONAL: make heatmap from an a priori list of genes ----
#read in a text file containing the genes (with expression data) you want to include in the heatmap
mySelectedGenes <- read_tsv("path/to/file/with/selected/genes/with/data")
#rather than reading a file in, this tibble could also come from dplyr operations in step 3 script
#convert to a matrix so you can carry out clustering
mySelectedGenes.matrix <- as.matrix(mySelectedGenes)
#you may (or may not) want to cluster your selected genes
hr <- hclust(as.dist(1-cor(t(mySelectedGenes.matrix), method="pearson")), method="complete") #cluster rows by pearson correlation
hc <- hclust(as.dist(1-cor(mySelectedGenes.matrix, method="spearman")), method="average") #cluster columns by spearman correlation

#make heatmap
heatmap.2(mySelectedGenes.matrix, 
          Rowv=NA, Colv=NA, 
          col=myheatcol, 
          scale="row", density.info="none", 
          trace="none", labCol=NA, 
          cexRow=1.5, cexCol=1, margins=c(8,20), key = F)


# the essentials ----
library(tidyverse)
library(gplots)
library(RColorBrewer)
myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") #cluster rows by pearson correlation
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20))

modulePick <- 2 
myModule_up <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub_up <- hclust(as.dist(1-cor(t(myModule_up), method="pearson")), method="complete") 

heatmap.2(myModule_up, 
          Rowv=as.dendrogram(hrsub_up), 
          Colv=NA, 
          labRow = NA,
          col=myheatcolors, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20))

modulePick <- 1 
myModule_down <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub_down <- hclust(as.dist(1-cor(t(myModule_down), method="pearson")), method="complete") 

heatmap.2(myModule_down, 
          Rowv=as.dendrogram(hrsub_down), 
          Colv=NA, 
          labRow = NA,
          col=myheatcolors, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20))

