# Load packages ------
library(dplyr)
library(readr) 
library(Biostrings) 
library(tximport)
library(ensembldb) 
library(EnsDb.Mmusculus.v79) 
library(ggplot2)
library(limma)
library(edgeR)
library(gplots) 
library(RColorBrewer)
library(limma)

targets <- read.table("studyDesign_Hackdash3.txt", row.names=NULL, header = T, as.is = T)
groups <- factor(paste(targets$genotype, targets$treatment, sep = "."))
sampleNames <- paste(targets$genotype, targets$treatment,targets$replicate, sep = ".")
design <- model.matrix(~0+groups)

Tx <- transcripts(EnsDb.Mmusculus.v79,columns=c(listColumns(EnsDb.Mmusculus.v79,"tx"), "gene_name"))
Tx <- as.data.frame(Tx)

Tx <- dplyr::rename(Tx, target_id = tx_id)
row.names(Tx) <- NULL
Tx <- Tx[,c(6,12)]

files <- file.path(targets$sample, "abundance.tsv")
all(file.exists(files))

Txi_gene <- tximport(files, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, 
                     countsFromAbundance = "no")

myTPM_raw <- Txi_gene$abundance
counts_filt <- Txi_gene$counts[rowSums(Txi_gene$counts)>=10,]

colnames(Txi_gene$counts) <- sampleNames
colnames(counts_filt) <- sampleNames

myTPM <- Txi_gene$abundance
distance <- dist(t(myTPM), method="maximum") 
clusters <- hclust(distance, method = "complete")
plot(clusters, labels=sampleNames, hang = -1)

# Pricipal component analysis (PCA) -------------
pca.res <- prcomp(t(myTPM), scale.=F, retx=T)
head(pca.res)
ls(pca.res)
summary(pca.res) 
head(pca.res$rotation)
head(pca.res$x)
plot(pca.res, las=1)
pc.var<-pca.res$sdev^2 
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per

# visualize your PCA result ------------------
#lets first plot any two PCs against each other
#turn your scores for each gene into a data frame
pcdata.frame <- as.data.frame(pca.res$x)
ggplot(pcdata.frame, aes(x=PC1, y=PC2, color=factor(groups))) +
  geom_point(size=5) +
  theme(legend.position="right")

# create a PCA 'small multiples' chart ----
# this is another way to view PCA to understand impact of each variable on each pricipal component
melted <- cbind(groups, melt(pca.res$x[,1:4]))
head(melted)
#look at your 'melted' data
ggplot(melted) +
  geom_bar(aes(x=Var1, y=value, fill=groups), stat="identity") +
  facet_wrap(~Var2)

# Normalization ------
myDGEList <- DGEList(counts_filt)
normData <- voom(myDGEList, design, plot = TRUE)
colnames(design) <- levels(groups)
# fit a linear model to your data
fit <- lmFit(normData, design)

# set up a contrast matrix based on the pairwise comparisons of interest ------
contrast.matrix.WTvsKO <- makeContrasts(unstim = WT.US-KO.US, FiveHr =WT.aCD328 - KO.aCD328,levels=design)
contrast.matrix.Tx <- makeContrasts(aCD328vsUS_WT = WT.aCD328 - WT.US,aCD328vsWT_KO = KO.aCD328 - KO.US,   levels=design)

# extract the linear model fit for the contrast matrix that you just defined above -----
fits.WTvsKO<- contrasts.fit(fit, contrast.matrix.WTvsKO)
fits.Tx<- contrasts.fit(fit, contrast.matrix.Tx)

#get bayesian stats for your linear model fit
ebFit.WTvsKO <- eBayes(fits.WTvsKO)
ebFit.Tx <- eBayes(fits.Tx)

# use decideTests to pull out the DEGs and make Venn Diagram ----
results.WTvsKO <- decideTests(ebFit.WTvsKO, method="global", adjust.method="BH", p.value=0.05, lfc=0.59)
results.Tx <- decideTests(ebFit.Tx, method="global", adjust.method="BH", p.value=0.01, lfc=1)

# take a look at what the results of decideTests looks like
vennDiagram(results.WTvsKO, include="both")
vennDiagram(results.Tx, include="both")

# retrieve expression data for your DEGs ----
colnames(normData$E) <- sampleNames

diffData.WTvsKO.unstim <- normData$E[results.WTvsKO[,1] != 0,]
diffData.WTvsKO.FiveHr <- normData$E[results.WTvsKO[,2] != 0,]
diffData.WT.aCD328_US <- normData$E[results.Tx[,1] != 0,]
diffData.KO.aCD328_US <- normData$E[results.Tx[,2] != 0,]

# Choose your color pallette ----
myheatcol <- greenred(75)

# generate a heatmap of DEGs ----
#begin by clustering the genes (rows) in each set of differentially expressed genes
hr.WTvsKO.unstim <- hclust(as.dist(1-cor(t(diffData.WTvsKO.unstim), method="pearson")), method="complete") #cluster rows by pearson correlation
hr.WTvsKO.FiveHr <- hclust(as.dist(1-cor(t(diffData.WTvsKO.FiveHr), method="pearson")), method="complete") #cluster rows by pearson correlation

mycl.WTvsKO.unstim <- cutree(hr.WTvsKO.unstim, k=2)
mycl.WTvsKO.FiveHr <- cutree(hr.WTvsKO.FiveHr, k=2)

#now assign a color to each cluster (makes it easy to identify and manipulate)
mycolhc.WTvsKO.unstim <- rainbow(length(unique(mycl.WTvsKO.unstim)), start=0.1, end=0.9) 
mycolhc.WTvsKO.unstim <- mycolhc.WTvsKO.unstim[as.vector(mycl.WTvsKO.unstim)] 

mycolhc.WTvsKO.FiveHr <- rainbow(length(unique(mycl.WTvsKO.FiveHr)), start=0.1, end=0.9) 
mycolhc.WTvsKO.FiveHr <- mycolhc.WTvsKO.FiveHr[as.vector(mycl.WTvsKO.FiveHr)] 

#plot the hclust results as a heatmap
heatmap.2(diffData.WTvsKO.unstim, Rowv=as.dendrogram(hr.WTvsKO.unstim), Colv=NA, 
          col=myheatcol, scale="row", labRow=NA,
          density.info="none", trace="none", RowSideColors=mycolhc.WTvsKO.unstim, 
          cexRow=1, cexCol=1, margins=c(10,10)) 

heatmap.2(diffData.WTvsKO.FiveHr, Rowv=as.dendrogram(hr.WTvsKO.FiveHr), Colv=NA, 
          col=myheatcol, scale="row", labRow=NA,
          density.info="none", trace="none", RowSideColors=mycolhc.WTvsKO.FiveHr, 
          cexRow=1, cexCol=1, margins=c(10,10)) 


# pull out clusters of co-regulated genes ----
# view your color assignments for the different clusters
names(mycolhc.WTvsKO.unstim) <- names(mycl.WTvsKO.unstim) 
barplot(rep(10, max(mycl.WTvsKO.unstim)),
        col=unique(mycolhc.WTvsKO.unstim[hr.WTvsKO.unstim$labels[hr.WTvsKO.unstim$order]]), 
        horiz=T, names=unique(mycl.WTvsKO.unstim[hr.WTvsKO.unstim$order]))

#choose a cluster(s) of interest by selecting the corresponding number based on the previous graph
clid <- c(1)
ysub <- diffData.WTvsKO.unstim[names(mycl.WTvsKO.unstim[mycl.WTvsKO.unstim%in%clid]),] 
hrsub <- hclust(as.dist(1-cor(t(ysub), method="pearson")), method="complete") 
clusterIDs <- data.frame(Labels=rev(hrsub$labels[hrsub$order]))
clusterIDs <- as.vector(t(clusterIDs))

# Create heatmap for chosen sub-cluster.
heatmap.2(ysub, Rowv=as.dendrogram(hrsub), Colv=NA, col=myheatcol, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=mycolhc.WTvsKO.unstim[mycl.WTvsKO.unstim%in%clid], margins=c(10,10)) 

clusterSymbols <- data.frame(Labels=rev(hrsub$labels[hrsub$order]))
clusterSymbols <- as.vector(t(clusterSymbols))
clusterData <- diffData.WTvsKO.unstim[clusterSymbols,]
write.csv(clusterData,"WTvsKO.unstim_cluster1.csv")
