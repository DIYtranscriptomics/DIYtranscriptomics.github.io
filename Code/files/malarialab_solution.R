# read data into R (step 1 script) ----
# this chunk contains the minimal essential code from this script. Simply uncomment the lines below and run the code.
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
targets <- read_tsv("studydesign.txt")# read in your study design
path <- file.path(targets$sample, "abundance.tsv") # set file paths to your mapped data
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     txOut = TRUE, #determines whether your data represented at transcript or gene level
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

# Filter and normalize (step 2 script) ----
library(edgeR)
library(matrixStats)
library(cowplot)

sampleLabels <- targets$sample
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = -1, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=2 #user defined
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = -1, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = -1, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)


# PCA and summary tables (step 3 script) ---- 
library(DT)
library(gt)
library(plotly)

group <- targets$group
#let's fix the group variable, since it may be ordered wrong (lexiographically), and since R doesn't like values that start with a number
targets$group <- c("timepoint_00hr",
                   "timepoint_00hr",
                   "timepoint_08hr",
                   "timepoint_08hr",
                   "timepoint_16hr", 
                   "timepoint_16hr",
                   "timepoint_24hr", 
                   "timepoint_24hr",
                   "timepoint_32hr", 
                   "timepoint_32hr",
                   "timepoint_40hr",
                   "timepoint_40hr",
                   "timepoint_48hr",
                   "timepoint_48hr")

group <- targets$group
group <- factor(group)


pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color = group) +
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
                    timepoint_0hr = (totalRNA_0hr_rep1 + totalRNA_0hr_rep2)/2, 
                    timepoint_8hr = (totalRNA_8hr_rep1 + totalRNA_8hr_rep2)/2, 
                    timepoint_16hr = (totalRNA_16hr_rep1 + totalRNA_16hr_rep2)/2, 
                    timepoint_24hr = (totalRNA_24hr_rep1 + totalRNA_24hr_rep2)/2, 
                    timepoint_32hr = (totalRNA_32hr_rep1 + totalRNA_32hr_rep2)/2, 
                    timepoint_48hr = (totalRNA_48hr_rep1 + totalRNA_48hr_rep2)/2, 
                    #now make columns comparing each of the averages above that you're interested in
                    LogFC = (timepoint_32hr - timepoint_0hr)) %>% #note that this is the first time you've seen the 'pipe' operator
  mutate_if(is.numeric, round, 2)

datatable(mydata.df[,c(1,16,20,22)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100")))

#sort the above table by LogFC. Top ranked genes should be PF3D7_0102200.1.  What does a quick lit/google search indicate about this gene?

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

# Different gene expression analysis (step 5 script) ----
group <- factor(targets$group)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

# Model mean-variance trend and fit linear model to data
# Use VOOM function from Limma package to model the mean-variance relationship
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)
# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design)

# Contrast matrix
contrast.matrix <- makeContrasts(DEGs_08hr = timepoint_08hr - timepoint_00hr,
                                 DEGs_16hr = timepoint_16hr - timepoint_00hr,
                                 DEGs_24hr = timepoint_24hr - timepoint_00hr,
                                 DEGs_32hr = timepoint_32hr - timepoint_00hr,
                                 DEGs_40hr = timepoint_40hr - timepoint_00hr,
                                 DEGs_48hr = timepoint_48hr - timepoint_00hr,
                                 levels=design)

# extract the linear model fit 
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)
#write.fit(ebFit, file="lmfit_results.txt")

# TopTable to view DEGs
myTopHits <- topTable(ebFit, adjust ="BH", coef=6, number=10, sort.by="logFC")

# convert to a tibble
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

gt(myTopHits.df)

# Volcano Plots
# rerun toptable with high number of genes selected
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=10000, sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")
# now plot
vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Malaria erythrocytic cycle",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# Now make the volcano plot above interactive with plotly
ggplotly(vplot)

# decideTests to pull out the DEGs and make Venn Diagram
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=1)

# take a look at what the results of decideTests looks like
head(results)
summary(results)
vennDiagram(results, include="up")

# retrieve expression data for your DEGs
head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0 | results[,3] !=0 | results[,4] !=0 | results[,5] !=0 | results[,6] !=0,]
head(diffGenes)
dim(diffGenes)
#convert your DEGs to a dataframe using as_tibble
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

# create interactive tables to display your DEGs
datatable(diffGenes.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: DEGs in cutaneous leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)

#write your DEGs to a file
write_tsv(diffGenes.df,"DiffGenes.txt") #NOTE: this .txt file can be directly used for input into other clustering or network analysis tools (e.g., String, Clust (https://github.com/BaselAbujamous/clust, etc.)

# heatmaps and module identification
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

# this heatmap is a good example of when clustering columns might be problematic.  Repeat heatmap without clustering columns.  Which do you prefer and why?
# Just-in-time gene expression: 
# originally described here: https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0000005
# more recent analysis of this phenomenon here: https://www.science.org/doi/10.1126/science.aba4357


