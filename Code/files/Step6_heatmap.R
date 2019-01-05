# introduction to this script -----------
#this script creates heatmaps from your differentially expressed genes or transcripts

# load packages -----
library(tidyverse)
library(gplots) 
library(RColorBrewer)
library(limma)
library(heatmaply) #for making interactive heatmaps using plotly
library(d3heatmap) #for making interactive heatmaps using D3

# choose color pallette ----
#Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]).
myheatcolors1 <- greenred(75)
# a color-blind friendly pallete
myheatcolors2 <- colorRampPalette(colors=c("yellow","white","blue"))(100)

# cluster DEGs ----
#begin by clustering the genes (rows) in each set of differentially expressed genes
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") #cluster rows by pearson correlation
# hierarchical clustering is a type of unsupervised clustering. Related methods include K-means, SOM, etc 
# unsupervised methods are blind to sample/group identity
# in contrast, supervised methods 'train' on a set of labeled data.  
# supervised clustering methods include random forest, and artificial neural networks

#now cluster your samples (columns)
#we may not acutally use this clustering result, but it's good to have just in case
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete") #cluster columns by spearman correlation
#note: we use Spearman, instead of Pearson, for clustering samples because it gives equal weight to highly vs lowly expressed transcripts or genes

# Cut the resulting tree and create color vector for clusters.  
#Vary the cut height to give more or fewer clusters, or you the 'k' argument to force n number of clusters
#we'll look at these clusters in more detail later
clust.assign <- cutree(clustRows, k=2)

#now assign a color to each cluster (makes it easy to identify and manipulate)
module.color <- rainbow(length(unique(clust.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(clust.assign)] 

# produce a static heatmap of DEGs ----
#plot the hclust results as a heatmap
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=NA,
          RowSideColors=module.color,
          col=myheatcolors2, scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20)) 

#what do the colors represent in this heatmap?
#what happens when you change scale=NULL

# make your heatmap interactive! ----
#first, we'll make an interactive heatmap using plotly (https://plot.ly/)
heatmaply(diffGenes,
          colors = myheatcolors2,
          Rowv=as.dendrogram(clustRows),
          RowSideColors=module.color,
          #showticklabels=c(FALSE,FALSE),
          scale='row')

# now let's try using D3 to create an html widget version of our heatmap
d3heatmap(diffGenes,
          colors = myheatcolors,
          Rowv=as.dendrogram(clustRows),
          row_side_colors = module.color,
          scale='row')


# you can annotate samples with any metadata available in your study design file
color.map <- function(groups) { if (groups=="uninfected") "#FF0000" else if (groups=="infected") "#33A12B" else "#0000FF"}
color.map <- unlist(lapply(groups, color.map))

heatmap.2(diffGenes, Rowv=as.dendrogram(clustRows), Colv=NA,
          col=myheatcolors, scale="row", labRow=NA,
          density.info="none", trace="none", 
          RowSideColors=module.color, ColSideColors = color.map,
          cexRow=1, cexCol=1, margins=c(8,20)) 

# simplify heatmap ----
#notice that the heatmap includes ALL the columns from your dataset
#to simplify, average biological replicates
#then rerun the heatmap the script above using diffData.AVG as input instead of diffData
colnames(diffGenes) <- groups

#now an old function from the limma package to average your replicates 
diffGenes.AVG <- avearrays(diffGenes)

##alternatively, decide exactly which columns you want to show, and modify the heatmap accordingly
#this is how it would look using base R
#diffGenes.subset <- diffGenes[,c(1,4,7)]
##now repeat heatmap only on these selected columns

# view clusters of co-regulated genes ----
# view your color assignments for the different clusters
names(module.color) <- names(clust.assign) 
barplot(rep(10, max(clust.assign)),
        col=unique(module.color[clustRows$labels[clustRows$order]]), 
        horiz=T, names=unique(clust.assign[clustRows$order]))

#choose a cluster(s) of interest by selecting the corresponding number based on the previous graph
clust.pick <- 2 #use c function to grab more than one cluster from the heatmap.  e.g., c(1,2)
mycluster <- diffGenes[names(clust.assign[clust.assign%in%clust.pick]),] 
hrsub <- hclust(as.dist(1-cor(t(mycluster), method="pearson")), method="complete") 
clusterIDs <- data.frame(Labels=rev(hrsub$labels[hrsub$order]))
clusterIDs <- as.vector(t(clust.assign))

# Create heatmap for chosen sub-cluster.
heatmap.2(mycluster, 
          Rowv=as.dendrogram(hrsub), 
          Colv=NA, 
          col=myheatcolors2, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[clust.assign%in%clust.pick], margins=c(8,20)) 

# print out clusters for downstream analysis ----
#prints out genes in the order you see them in the cluster
clusterSymbols <- data.frame(Labels=rev(hrsub$labels[hrsub$order]))
clusterSymbols <- as.vector(t(clusterSymbols))
clusterData <- diffGenes[clusterSymbols,]
clusterData.df <- as_tibble(clusterData, rownames = "geneSymbol")
write_csv(clusterData.df,"Cluster_downRegulated.csv")

# OPTIONAL: make heatmap from an a priori list of genes ----
#read in a text file containing the genes (with expression data) you want to include in the heatmap
mySelected <- read_tsv("mySelectedGenes.txt")

mySelected.matrix <- as.matrix(mySelected)
#alternatively, you can go straight from a dataframe filtered by dplyr
myTPM.filter
rownames(myTPM.filter) <- myTPM.filter[,1]
myTPM.filter <- myTPM.filter[,-1]
myTPM.filter <- as.matrix(myTPM.filter)
#you may (or may not) want to cluster your selected genes
hr <- hclust(as.dist(1-cor(t(myTPM.filter), method="pearson")), method="complete") #cluster rows by pearson correlation
hc <- hclust(as.dist(1-cor(myTPM.filter, method="spearman")), method="average") #cluster columns by spearman correlation

#make heatmap
heatmap.2(myTPM.filter, Rowv=NA, Colv=NA, 
          col=myheatcol, scale="row", density.info="none", 
          trace="none", labCol=NA, 
          cexRow=1.5, cexCol=1, margins=c(8,20), key = F) # Creates heatmap for entire data set where the obtained clusters are indicated in the color bar.

