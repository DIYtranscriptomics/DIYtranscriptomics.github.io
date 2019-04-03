library(tidyverse)
library(gt)
library(hrbrthemes)

# read in data and study design ----
#if the file is in the same directory as our analysis script, then no path is necessary
data <- read_tsv("Schisto_Log2CPM.unfiltered.txt")

#our data is in the form of a dataframe, but we need a data matrix for PCA
#first, use subsetting operation to get rid of columns that are non numeric
#then use the 'as.matrix' function to create a new data matrix from this numeric subset of the data
data.matrix <- data[,-c(1:3)]
data.matrix <- as.matrix(data.matrix)
studyDesign <- read_tsv("Schisto_studyDesign.txt")

# PCA ----
pca.res <- prcomp(t(data.matrix), scale.=F, retx=T)
summary(pca.res) # Prints variance summary for all principal components.
x <- pca.res$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
pca.res$x #$x shows you how much each sample influenced each PC (called 'loadings')
#note that these loadings have a magnitude and a direction (this is the basis for making a PCA plot)
pc.var<-pca.res$sdev^2 #sdev^2 gives you the eigenvalues
pc.per<-round(pc.var/sum(pc.var)*100, 1)

# Visualize PCA result ------------------
#lets first plot any two PCs against each other
#We know how much each sample contributes to each PC (loadings), so let's plot
pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df, aes(x=PC1, y=PC2, color=studyDesign$parasite_drugSensitivity_praziquantel)) +
  geom_point(size=3) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       subtitle="PC1 reflects parasite sex",
       caption=paste0("produced on ", Sys.time())) +
  theme_ipsum_rc()

#now recolor according to another variables and PCs
pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df, aes(x=PC3, y=PC4, color=studyDesign$parasiteStrain)) +
  geom_point(size=3) +
  xlab(paste0("PC3 (",pc.per[3],"%",")")) + 
  ylab(paste0("PC4 (",pc.per[4],"%",")")) +
  labs(title="PCA plot",
       subtitle="PC3/4 separate samples based on parasite strain",
       caption=paste0("produced on ", Sys.time())) +
  theme_ipsum_rc()

# dplyr to ID top 10 genes ----
# begin by selecting the columns based on sex, strain and timepoint
# lots of ways to do this, but an easy one is to use dplyr select and the 'starts_with' function
data.subset <- data %>%
  dplyr::select(starts_with("F24h_LE_"), starts_with("FCtl_LE_")) %>%
  dplyr::mutate(geneID = data$wbps_transcript_id)

#now use dplyr mutate to create averages for control and 24hr, as well as a logFC column
data.subset <- data.subset %>%
  dplyr::mutate(F24h_LE_AVG = (F24h_LE_1 + F24h_LE_2 + F24h_LE_3)/3,  
                FCtl_LE_AVG = (FCtl_LE_1 + FCtl_LE_2 + FCtl_LE_3)/3, 
                LogFC_24hr.vs.Ctl = F24h_LE_AVG - FCtl_LE_AVG)

#now use dplyr filter and arrange to pick and rank the top 10 genes based on logFC
data.subset <- data.subset %>%
  dplyr::arrange(desc(LogFC_24hr.vs.Ctl)) %>%
  dplyr::top_n(10)


#now cleaning up a bit
data.subset <- data.subset %>%
  dplyr::select(geneID, FCtl_LE_AVG, F24h_LE_AVG, LogFC_24hr.vs.Ctl)

#make the table
gt(data.subset)
