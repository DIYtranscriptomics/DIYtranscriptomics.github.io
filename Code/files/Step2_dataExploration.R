# Introduction to this script ----
# now that you've read your transcript-level or gene-level data into R, you're ready to begin exploring the structure of the data
# goal of this script is to using multivariate statisical approaches to explore the structure of your data
# recall that your abundance data are TPM, while the counts are read counts mapping to each gene or transcript

# Load packages -----
library(ggplot2)
library(reshape2)
library(trelliscopejs)
library(genefilter)
library(limma)
library(edgeR)
library(Biobase)

# Read in your study design ----
targets <- read.table("Crypto_studyDesign.txt", row.names=NULL, header = T, as.is = T)
groups <- targets$treatment
groups <- factor(groups)
sampleLabels <- targets$sample
#what would you do if there more variables than just treatment in your experiment?

# Examine your data up to this point ----
myTPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts
# Do you think counts or TPM are appropriate for multivariate statistical analysis....why?
# let's quickly graph both matrices to see what we're dealing with

# Take a look at the heteroskedasticity of the data ----
# first, calculate row means and standard deviations for each transcript or gene and add these to your data frame
myTPM.stats <- transform(myTPM, SD=rowSds(myTPM), AVG=rowMeans(myTPM), MED=rowMedians(myTPM))
myCounts.stats <- transform(myCounts, SD=rowSds(myCounts), AVG=rowMeans(myCounts), MED=rowMedians(myCounts))

head(counts.stats)
ggplot(TPM.stats, aes(x=SD, y=MED)) +
  geom_point(shape=1) +
  geom_point(size=4)
# how might you expect that counts and TPM compare if used as input for PCA analysis?
# how would these graphs change if you log2 converted the data?

# Make a DGElist from your counts ----
DGEList <- DGEList(Txi_gene$counts)
# use the 'cpm' function from EdgeR to get counts per million
cpm <- cpm(DGEList) 
log2.cpm <- cpm(DGEList, log=TRUE)

# Take a look at the distribution of the Log2 CPM
library(RColorBrewer)
nsamples <- ncol(log2.cpm)
myColors <- brewer.pal(nsamples, "Paired")
boxplot(log2.cpm, las=2, cex=1, col = myColors, names = sampleLabels, main="non-normalized log2 cpm")
#what do you see as some potential issues here?
# Filter your data ----
#first, take a look at how many genes or transcripts have no read counts at all
table(rowSums(DGEList.norm$counts==0)==9)
# now set some cut-off to get rid of genes/transcripts with low counts
keepers <- rowSums(cpm>1)>=3
DGEList.filtered <- DGEList[keepers,]
dim(DGEList.filtered)

log2.cpm.filtered <- cpm(DGEList.filtered, log=TRUE)
boxplot(log2.cpm.filtered, 
        las=2, cex=1, 
        col = myColors, 
        names = sampleLabels, 
        main="non-normalized log2 cpm")

# Normalize your data ----
DGEList.filtered.norm <- calcNormFactors(DGEList.filtered, method = "TMM")

# use the 'cpm' function from EdgeR to get counts per million from you 
log2.cpm.filtered.norm <- cpm(DGEList.filtered.norm, log=TRUE)
boxplot(log2.cpm.filtered.norm, 
        las=2, cex=1, 
        col = myColors, 
        names = sampleLabels, 
        main="normalized log2 cpm")


# Hierarchical clustering ---------------
#hierarchical clustering can only work on a data matrix, not a data frame
#try using filtered and unfiltered data...how does this change the results
distance <- dist(t(log2.cpm.filtered.norm), method="maximum") #other dist methods are "maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "complete") #other methods are ward.D, ward.D2, single, complete, average
plot(clusters, labels=sampleLabels)

# Pricipal component analysis (PCA) -------------
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
#look at pca.res in environment
class(pca.res)
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
head(pca.res$rotation) #$rotation shows you how much each TRANSCRIPT influenced each PC (called 'loadings')
head(pca.res$x) #$x shows you how much each SAMPLE influenced each PC (called 'scores')
plot(pca.res, las=1)
pc.var<-pca.res$sdev^2 #sdev^2 gives you the eigenvalues
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per

# Visualize your PCA result ------------------
#lets first plot any two PCs against each other
#turn your scores for each gene into a data frame
data.frame <- as.data.frame(pca.res$x)
ggplot(data.frame, aes(x=PC1, y=PC2, color=groups)) +
  geom_point(size=5) +
  theme(legend.position="right")

# Create a PCA 'small multiples' chart ----
# this is another way to view PCA to understand impact of each variable on each pricipal component
melted <- cbind(groups, melt(pca.res$x[,1:4]))
head(melted)
#look at your 'melted' data
ggplot(melted) +
  geom_bar(aes(x=Var1, y=value, fill=groups), stat="identity") +
  facet_wrap(~Var2)
#facet_trelliscope(~Var2)
# Explore other plotting options on your own ----
# T-sne, NMDS, reactive graphs, 3D interactives, etc
library(scatterD3) #makes nice interactive 2D plots using the D3 engine
library(rgl)
library(pca3d) # good for 3D PCA plots
library(RColorBrewer)
test <- cbind(pca.res$x, c(1,1,1,2,2,2,3,3,3))
colnames(test)[10] <- "colvector"
scatterD3(test, 
          x = PC1, y = PC2,
          col_var = colvector)
pca3d(pca.res$x[,1:3], 
      radius = 3,
      col = c("blue", "blue", "blue", "red", "red", "red", "green", "green", "green"))
