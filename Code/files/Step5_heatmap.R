# Introduction to this script -----------
#this script creates heatmaps from your differentially expressed genes or transcripts

# Load packages -----
library(gplots) 
library(RColorBrewer)
library(limma)

# Choose your color pallette ----
#Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]).
myheatcol <- greenred(75)
# a color-blind friendly pallete
myheatcol <- colorRampPalette(colors=c("yellow","white","blue"))(100)

# generate a heatmap of DEGs ----
#begin by clustering the genes (rows) in each set of differentially expressed genes
hr <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") #cluster rows by pearson correlation

#now cluster your samples (columns)
#we may not acutally use this clustering result, but it's good to have just in case
hc <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete") #cluster columns by spearman correlation

# Cut the resulting tree and create color vector for clusters.  
#Vary the cut height to give more or fewer clusters, or you the 'k' argument to force n number of clusters
#we'll look at these clusters in more detail later
mycl <- cutree(hr, k=2)

#now assign a color to each cluster (makes it easy to identify and manipulate)
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9) 
mycolhc <- mycolhc[as.vector(mycl)] 

#plot the hclust results as a heatmap
# first, a heatmap of genes regulated by LPS in a WT background
heatmap.2(diffGenes, Rowv=as.dendrogram(hr), Colv=NA, 
          col=myheatcol, scale="row", labRow=NA,
          density.info="none", trace="none", RowSideColors=mycolhc, 
          cexRow=1, cexCol=1, margins=c(8,20)) 

# edit heatmap to simplify----
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

# pull out clusters of co-regulated genes ----
# view your color assignments for the different clusters
names(mycolhc) <- names(mycl) 
barplot(rep(10, max(mycl)),
        col=unique(mycolhc[hr$labels[hr$order]]), 
        horiz=T, names=unique(mycl[hr$order]))

#choose a cluster(s) of interest by selecting the corresponding number based on the previous graph
#first cluster to pick are the genes in KO cells that lose LPS-inducibility
clid <- c(2)
ysub <- diffGenes[names(mycl[mycl%in%clid]),] 
hrsub <- hclust(as.dist(1-cor(t(ysub), method="pearson")), method="complete") 
clusterIDs <- data.frame(Labels=rev(hrsub$labels[hrsub$order]))
clusterIDs <- as.vector(t(clusterIDs))

# Create heatmap for chosen sub-cluster.
heatmap.2(ysub, Rowv=as.dendrogram(hrsub), Colv=NA, col=myheatcol, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=mycolhc[mycl%in%clid], margins=c(8,20)) 

# print out clusters for downstream analysis ----
#prints out genes in the order you see them in the cluster
clusterSymbols <- data.frame(Labels=rev(hrsub$labels[hrsub$order]))
clusterSymbols <- as.vector(t(clusterSymbols))
clusterData <- diffGenes[clusterSymbols,]
write.table(clusterData,"Cluster1.xls", sep="\t", quote=FALSE)

# OPTIONAL: make heatmap from an a priori list of genes ----
#read in a text file containing the genes (with expression data) you want to include in the heatmap
mySelected <- read.delim("mySelectedGenes.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, row.names=1)
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

