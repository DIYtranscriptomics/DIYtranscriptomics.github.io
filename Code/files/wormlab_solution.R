# SET UP - read in data, study design, and capture variables of interest ----
library(tidyverse)
library(tximport)
targets <- read_tsv("studyDesign.txt")
sampleLabels <- targets$sample
sex <- factor(targets$sex)
strain <- factor(targets$strain)
timepoint <- factor(targets$timpoint)
treatment <- factor(targets$treatment)
drugSensitivity <- factor(targets$drugSensitivity)

# read in Kallisto data
path <- file.path(targets$sample, "abundance.h5")
Txi_gene <- tximport(path,
                     type = "kallisto",
                     txOut = TRUE,
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)



# FILTER, NORMALIZE, PLOT (question 1)----
library(reshape2)
library(edgeR)
library(matrixStats)
library(cowplot)
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = 2:145, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression")

ggplot(log2.cpm.df.pivot, aes(x=samples, y=expression, fill=samples)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.y = "median",
               geom = "point",
               shape = 124,
               size = 6,
               color = "black",
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time()))

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=3 #user defined
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)

log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = 2:145, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

ggplot(log2.cpm.filtered.df.pivot, aes(x=samples, y=expression, fill=samples)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.y = "median",
               geom = "point",
               shape = 124,
               size = 6,
               color = "black",
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time()))

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
# pivot this NORMALIZED data, just as you did earlier
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = 2:145, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

ggplot(log2.cpm.filtered.norm.df.pivot, aes(x=samples, y=expression, fill=samples)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.y = "median",
               geom = "point",
               shape = 124,
               size = 6,
               color = "black",
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time()))


# EXPLORE BY PCA (question 2) ----
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var<-pca.res$sdev^2
pc.per<-round(pc.var/sum(pc.var)*100, 1)

pca.res.df <- as_tibble(pca.res$x)

#Here, I chose to mape point color to the worm strain, and point shape to
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color=sex) +
  geom_point(size=4) +
  # geom_label() +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) +
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  # coord_fixed() +
  theme_bw()

pca.res.df <- pca.res$x[,1:6] %>%
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = targets$sex)

pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC6, # column names to be stored as a SINGLE variable
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

# interpretation of PCA
# PC1 is sex and explains nearly 50% of the variance in the total dataset
# PC2 and PC4 are strain and collectively explain about 20% of the variance
# PC3 seems to be linked to timepoint

# let's confirm our interpretations from the small multiples by returning to our original PCA
#Let's look at PC2 vs PC4
ggplot(pca.res.df) +
  aes(x=PC2, y=PC4, color=strain) +
  geom_point(size=4) +
  # geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[2],"%",")")) +
  ylab(paste0("PC2 (",pc.per[4],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  # coord_fixed() +
  theme_bw()


# DATA GYMNASTICS WITH DPLYR (question 3) ----
# begin by selecting the columns based on sex, strain and timepoint
# lots of ways to do this, but an easy one is to use dplyr select and the 'starts_with' function
data.subset <- log2.cpm.filtered.norm.df %>%
  dplyr::select(geneID, starts_with("F24h_LE_"), starts_with("FCtl_LE_"))

#now use dplyr mutate to create averages for control and 24hr, as well as a logFC column
data.subset <- data.subset %>%
  dplyr::mutate(F24h_LE_AVG = (F24h_LE_1 + F24h_LE_2 + F24h_LE_3)/3,
                FCtl_LE_AVG = (FCtl_LE_1 + FCtl_LE_2 + FCtl_LE_3)/3,
                LogFC_24hr.vs.Ctl = F24h_LE_AVG - FCtl_LE_AVG)

#now use dplyr filter and arrange to pick and rank the top 10 genes based on logFC
data.subset <- data.subset %>%
  dplyr::mutate_if(is.numeric, round, 2) %>%
  dplyr::arrange(desc(LogFC_24hr.vs.Ctl)) %>%
  dplyr::top_n(10) %>%
  dplyr::select(geneID, FCtl_LE_AVG, F24h_LE_AVG, LogFC_24hr.vs.Ctl)

# it's easy to make nice tables with the gt package
library(gt)
gt(data.subset)

# ANY INTERSTING BIOLOGY IN THIS LIST OF GENES (Bonus question) ----
# multiple genes in the top 10 list (e.g. Smp_138080.1) belong to a unique class of genes called Micro Exon Genes (MEGs)
# MEGs have many exons..usually very small, and always of a length that is divisble by 3.
# as a result, the transcriptional machinery easily 'slips' past an exon, but the transcript remains in-frame
# the result is genetic diversity without need to rely on alternative splicing
# MEGs also happen to be secreted, which may be one reason the parasite wants to easily and quickly modify the protein
# super interesting that these genes seem to be induced by drug treatment.
# could this be one way the parasite tries to evade the immune response during treatment?
