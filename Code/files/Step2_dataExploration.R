# Introduction to this script ----
# now that you've read your transcript-level or gene-level data into R, you're ready to begin exploring the structure of the data
# goal of this script is to using multivariate statisical approaches to explore the structure of your data
# begin by taking a look at the 'abundance' element in the 'Txi' or 'Txi_gene' object you created at the end of the last class
# this is your Transcripts per Million (TPM) measurements
head(Txi_gene$abundance)
myTPM <- Txi_gene$abundance
head(myTPM)
# Load packages -----
library(ggplot2)
library(reshape2)
library(trelliscopejs)

# set-up study design ----
targets <- read.table("Crypto_studyDesign.txt", row.names=NULL, header = T, as.is = T)
groups <- targets$treatment
groups <- factor(groups)
sampleLabels <- targets$sample
#what would you do if there more variables than just treatment in your experiment?

# hierarchical clustering ---------------
#hierarchical clustering can only work on a data matrix, not a data frame
#try using filtered and unfiltered data...how does this change the results
distance <- dist(t(myTPM), method="maximum") #other dist methods are "maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "complete") #other methods are ward.D, ward.D2, single, complete, average
plot(clusters, labels=sampleLabels)

# Pricipal component analysis (PCA) -------------
pca.res <- prcomp(t(myTPM), scale.=F, retx=T)
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

# visualize your PCA result ------------------
#lets first plot any two PCs against each other
#turn your scores for each gene into a data frame
data.frame <- as.data.frame(pca.res$x)
ggplot(data.frame, aes(x=PC1, y=PC3, color=groups)) +
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
  #facet_trelliscope(~Var2)
# explore other plotting options on your own ----
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