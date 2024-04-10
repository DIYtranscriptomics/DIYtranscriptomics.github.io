# Introduction to this script ----
# for the purposes of this script we'll want several data objects generated in previous scripts, including:
# 1) your normalized filtered expression data, in the form of a data matrix with symbols as rownames.
# 2) your study design file
# 3) your contrast matrix that lays out the pairwise comparisons you're interested in testing
# 4) Individual signatures or 'collections' of signatures to test for enrichment in your data.
# These signatures can be downloaded from gene signature databases such as MSigDB
# Signatures can also be custom made based on your interests.
# Signatures can also be pulled from R/Bioconductor as described below

# Load packages ----
library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots

# Carry out GO enrichment using gProfiler2 ----
# use topTable result to pick the top genes for carrying out a Gene Ontology (GO) enrichment analysis
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=50, sort.by="logFC")
# use the 'gost' function from the gprofiler2 package to run GO enrichment analysis
gost.res <- gost(rownames(myTopHits), organism = "hsapiens", correction_method = "fdr")
# produce an interactive manhattan plot of enriched GO terms
gostplot(gost.res, interactive = TRUE, capped = TRUE)
#set interactive=FALSE to get plot for publications
mygostplot <- gostplot(gost.res, interactive = FALSE, capped = TRUE)
# produce a publication quality static manhattan plot with specific GO terms highlighted
publish_gostplot(
  mygostplot, #your static gostplot from above
  highlight_terms = c("REAC:R-HSA-9662851", "REAC:R-HSA-9824443", "REAC:R-HSA-9664407", "REAC:R-HSA-9658195", "REAC:R-HSA-9664417"),
  filename = NULL,
  width = NA,
  height = NA)

# now repeat the above steps using only genes from a single module from the step 6 script, by using `rownames(myModule)`
# what is value in breaking up DEGs into modules for functional enrichment analysis?

# Perform GSEA using clusterProfiler ----
# there are a few ways to get msigDB collections into R
# option1: download directly from msigdb and load from your computer
# can use the 'read.gmt' function from clusterProfiler package to create a dataframe,
# alternatively, you can read in using 'getGmt' function from GSEABase package if you need a GeneSetCollection object
C2CP <- read.gmt("/Users/danielbeiting/Dropbox/MSigDB/c2.cp.v2023.2.Hs.symbols.gmt")

# option2: use the msigdb package to access up-to-date collections
# this option has the additional advantage of providing access to species-specific collections
# are also retrieved as tibbles
msigdbr_species()
hs_gsea <- msigdbr(species = "Homo sapiens") #gets all collections/signatures with human gene IDs
#take a look at the categories and subcategories of signatures available to you
hs_gsea %>%
  dplyr::distinct(gs_cat, gs_subcat) %>%
  dplyr::arrange(gs_cat, gs_subcat)

# choose a specific msigdb collection/subcollection
# since msigdbr returns a tibble, we'll use dplyr to do a bit of wrangling
hs_gsea_c2 <- msigdbr(species = "Homo sapiens", # change depending on species your data came from
                      category = "C2") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature

# Now that you have your msigdb collections ready, prepare your data
# grab the dataframe you made in step3 script
# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)
# construct a named vector
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
set.seed(123) #set a random seed so that we can reproducible ordering for our GSEA results below
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=C2CP, verbose=FALSE) #could replace C2CP with hs_gsea_c2 object you retrieved from msigdb above
myGSEA.df <- as_tibble(myGSEA.res@result)

# view results as an interactive table
datatable(myGSEA.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:10), digits=2)
# create enrichment plots using the enrichplot package
gseaplot2(myGSEA.res,
          geneSetID = c(9, 23, 28, 35), #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          #title = myGSEA.res$Description[14]
          ) #can also turn off this title

# add a variable to this result that matches enrichment direction with phenotype
myGSEA.df <- myGSEA.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "disease",
    NES < 0 ~ "healthy"))

# create 'bubble plot' to summarize y signatures across x phenotypes
ggplot(myGSEA.df[1:20,], aes(x=phenotype, y=ID)) +
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()

# Competitive GSEA using CAMERA----
# for competitive tests the null hypothesis is that genes in the set are, at most, as often differentially expressed as genes outside the set
# first let's create a few signatures to test in our enrichment analysis
mySig <- rownames(myTopHits) #if your own data is from mice, wrap this in 'toupper()' to make gene symbols all caps
mySig2 <- sample((rownames(v.DEGList.filtered.norm$E)), size = 50, replace = FALSE)
collection <- list(real = mySig, fake = mySig2)
# now test for enrichment using CAMERA
camera.res <- camera(v.DEGList.filtered.norm$E, collection, design, contrast.matrix[,1])
camera.df <- as_tibble(camera.res, rownames = "setName")
camera.df

# Let's compare this with what happens when we use ROAST (also from the Limma package)
# ROAST uses a self-contained approach, which evaluates the null hypothesis that no genes in the set are differentially expressed
mroast(v.DEGList.filtered.norm$E, collection, design, contrast=1) #mroast adjusts for multiple testing

# now repeat CAMERA with an actual gene set collection
# camera requires collections to be presented as a list, rather than a tibble, so we must read in our signatures using the 'getGmt' function
C2CP <- getGmt("/Users/danielbeiting/Dropbox/MSigDB/c2.cp.v2023.2.Hs.symbols.gmt", geneIdType=SymbolIdentifier())
#extract as a list
C2CP <- geneIds(C2CP)
camera.res <- camera(v.DEGList.filtered.norm$E, C2CP, design, contrast.matrix[,1])
camera.df <- as_tibble(camera.res, rownames = "setName")
camera.df

# filter based on FDR and display as interactive table
camera.df <- dplyr::filter(camera.df, FDR<=0.01)

datatable(camera.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2,4,5), digits=2)

#as before, add a variable that maps up/down regulated pathways with phenotype
camera.df <- camera.df %>%
  mutate(phenotype = case_when(
    Direction == "Up" ~ "disease",
    Direction == "Down" ~ "healthy"))

#easy to filter this list based on names of signatures using 'str_detect'
#here is an example of filtering to return anything that has 'CD8' or 'CYTOTOX' in the name of the signature
camera.df.sub <- camera.df %>%
  dplyr::filter(str_detect(setName, "CD8|CYTOTOX"))

# graph camera results as bubble chart
ggplot(camera.df[1:25,], aes(x=phenotype, y=setName)) +
  geom_point(aes(size=NGenes, color = Direction, alpha=-log10(FDR))) +
  theme_bw()

# Single sample GSEA using the GSVA package----
# the GSVA package offers a different way of approaching functional enrichment analysis.
# A few comments about the approach:
# In contrast to most GSE methods, GSVA performs a change in coordinate systems,
# transforming the data from a gene by sample matrix to a gene set (signature) by sample matrix.
# this allows for the evaluation of pathway enrichment for each sample.
# the method is both non-parametric and unsupervised
# bypasses the conventional approach of explicitly modeling phenotypes within enrichment scoring algorithms.
# focus is therefore placed on the RELATIVE enrichment of pathways across the sample space rather than the absolute enrichment with respect to a phenotype.
# however, with data with a moderate to small sample size (< 30), other GSE methods that explicitly include the phenotype in their model are more likely to provide greater statistical power to detect functional enrichment.

# be aware that if you choose a large MsigDB file here, this step may take a while
# we'll use the canonical pathways portion of C2 from MSigDB that we prepared earlier
# first, we build a GSVA parameter object
gsvapar <- gsvaParam(exprData = v.DEGList.filtered.norm$E,
                     geneSets = C2CP,
                     minSize=5, maxSize=500,
                     maxDiff=TRUE)
# now we run GSVA using this parameter
GSVA.res.C2CP <- gsva(gsvapar)
# Apply linear model to GSVA result
# now using Limma to find significantly enriched gene sets in the same way you did to find diffGenes
# this means you'll be using topTable, decideTests, etc
# note that you need to reference your design and contrast matrix here
fit.C2CP <- lmFit(GSVA.res.C2CP, design)
ebFit.C2CP <- eBayes(fit.C2CP)

# use topTable and decideTests functions to identify the differentially enriched gene sets
topPaths.C2CP <- topTable(ebFit.C2CP, adjust ="BH", coef=1, number=50, sort.by="logFC")
res.C2CP <- decideTests(ebFit.C2CP, method="global", adjust.method="BH", p.value=0.05, lfc=0.5)
# the summary of the decideTests result shows how many sets were enriched in induced and repressed genes in all sample types
summary(res.C2CP)

# pull out the gene sets that are differentially enriched between groups
diffSets.C2CP <- GSVA.res.C2CP[res.C2CP[,1] !=0,]
head(diffSets.C2CP)
dim(diffSets.C2CP)

# make a heatmap of differentially enriched gene sets
hr.C2CP <- hclust(as.dist(1-cor(t(diffSets.C2CP), method="pearson")), method="complete") #cluster rows by pearson correlation
hc.C2CP <- hclust(as.dist(1-cor(diffSets.C2CP, method="spearman")), method="complete") #cluster columns by spearman correlation

# Cut the resulting tree and create color vector for clusters.  Vary the cut height to give more or fewer clusters, or you the 'k' argument to force n number of clusters
mycl.C2CP <- cutree(hr.C2CP, k=2)
mycolhc.C2CP <- rainbow(length(unique(mycl.C2CP)), start=0.1, end=0.9)
mycolhc.C2CP <- mycolhc.C2CP[as.vector(mycl.C2CP)]

# assign your favorite heatmap color scheme. Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]). Type demo.col(20) to see more color schemes.
myheatcol <- colorRampPalette(colors=c("yellow","white","blue"))(100)
# plot the hclust results as a heatmap
heatmap.2(diffSets.C2CP,
          Rowv=as.dendrogram(hr.C2CP),
          Colv=NA,
          col=myheatcol, scale="row",
          density.info="none", trace="none",
          cexRow=0.9, cexCol=1, margins=c(10,14)) # Creates heatmap for entire data set where the obtained clusters are indicated in the color bar.
# just as we did for genes, we can also make an interactive heatmap for pathways
# you can edit what is shown in this heatmap, just as you did for your gene level heatmap earlier in the course
# the essentials (gProfiler2 for GO enrichment)----
# this code chunk assumes you already have a contrast matrix and ebFit object from the Step5_diffGenes script
# load packages
library(tidyverse)
library(limma)
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources

# use topTable result to pick the top genes for carrying out a Gene Ontology (GO) enrichment analysis
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=50, sort.by="logFC")
# use the 'gost' function from the gprofiler2 package to run GO enrichment analysis
gost.res <- gost(rownames(myTopHits), organism = "hsapiens", correction_method = "fdr")
# produce an interactive manhattan plot of enriched GO terms
gostplot(gost.res, interactive = TRUE, capped = TRUE)

# the essentials (clusterProfiler for GSEA)----
# this code chunk assumes you have a local copy of an msigdb database and a dataframe with logFCs from the Step3_multivariate script
# load packages
library(tidyverse)
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(enrichplot) # great for making the standard GSEA enrichment plots

# read in your signature database
C2CP <- read.gmt("/Users/danielbeiting/Dropbox/MSigDB/c2.cp.v2023.2.Hs.symbols.gmt")
# grab the dataframe you made in step3 script
# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)
# construct a named vector
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
set.seed(123) #set a random seed so that we can reproducible ordering for our GSEA results below
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=C2CP, verbose=FALSE)
myGSEA.df <- as_tibble(myGSEA.res@result)

# view results as an interactive table
datatable(myGSEA.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:10), digits=2)
# create enrichment plots using the enrichplot package
gseaplot2(myGSEA.res,
          geneSetID = c(14, 27, 28, 30), #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          #title = myGSEA.res$Description[14]
          ) #can also turn off this title

# add a variable to this result that matches enrichment direction with phenotype
myGSEA.df <- myGSEA.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "disease",
    NES < 0 ~ "healthy"))

# create 'bubble plot' to summarize y signatures across x phenotypes
ggplot(myGSEA.df[1:20,], aes(x=phenotype, y=ID)) +
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()

# the essentials (CAMERA for GSEA)----
# this chunk assumes you have a contrast matrix and design set-up (see step 5 script), as well as a signature database (I'll use my local one again)
# load packages
library(tidyverse)
library(limma)
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase

# read in signature database
C2CP <- getGmt("/Users/danielbeiting/Dropbox/MSigDB/c2.all.v2023.2.Hs.symbols.gmt", geneIdType=SymbolIdentifier())
C2CP <- geneIds(C2CP) #extract as a list
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)

camera.res <- camera(v.DEGList.filtered.norm$E, C2CP, design, contrast.matrix[,1])
camera.df <- as_tibble(camera.res, rownames = "setName")
camera.df

# filter based on FDR and display as interactive table
camera.df <- dplyr::filter(camera.df, FDR<=0.01)

datatable(camera.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2,4,5), digits=2)

#as before, add a variable that maps up/down regulated pathways with phenotype
camera.df <- camera.df %>%
  mutate(phenotype = case_when(
    Direction == "Up" ~ "disease",
    Direction == "Down" ~ "healthy"))

# graph camera results as bubble chart
ggplot(camera.df[1:25,], aes(x=phenotype, y=setName)) +
  geom_point(aes(size=NGenes, color = Direction, alpha=-log10(FDR))) +
  theme_bw()


# the essentials (GSVA) ----
# this chunk assumes you have a contrast matrix and design set-up (see step 5 script)
library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.

# read in your signature database
C2CP <- read.gmt("/Users/danielbeiting/Dropbox/MSigDB/c2.cp.v2023.2.Hs.symbols.gmt")

# first, we build a GSVA parameter object
gsvapar <- gsvaParam(exprData = v.DEGList.filtered.norm$E,
                     geneSets = C2CP,
                     minSize=5, maxSize=500,
                     maxDiff=TRUE)
# now we run GSVA using this parameter
GSVA.res.C2CP <- gsva(gsvapar)
# Apply linear model to GSVA result
# now using Limma to find significantly enriched gene sets in the same way you did to find diffGenes
# this means you'll be using topTable, decideTests, etc
# note that you need to reference your design and contrast matrix here
fit.C2CP <- lmFit(GSVA.res.C2CP, design)
ebFit.C2CP <- eBayes(fit.C2CP)

# use topTable and decideTests functions to identify the differentially enriched gene sets
topPaths.C2CP <- topTable(ebFit.C2CP, adjust ="BH", coef=1, number=50, sort.by="logFC")
res.C2CP <- decideTests(ebFit.C2CP, method="global", adjust.method="BH", p.value=0.05, lfc=0.5)
# the summary of the decideTests result shows how many sets were enriched in induced and repressed genes in all sample types
summary(res.C2CP)

# pull out the gene sets that are differentially enriched between groups
diffSets.C2CP <- GSVA.res.C2CP[res.C2CP[,1] !=0,]

# make a heatmap of differentially enriched gene sets
hr.C2CP <- hclust(as.dist(1-cor(t(diffSets.C2CP), method="pearson")), method="complete") #cluster rows by pearson correlation
hc.C2CP <- hclust(as.dist(1-cor(diffSets.C2CP, method="spearman")), method="complete") #cluster columns by spearman correlation

# Cut the resulting tree and create color vector for clusters.  Vary the cut height to give more or fewer clusters, or you the 'k' argument to force n number of clusters
mycl.C2CP <- cutree(hr.C2CP, k=2)
mycolhc.C2CP <- rainbow(length(unique(mycl.C2CP)), start=0.1, end=0.9)
mycolhc.C2CP <- mycolhc.C2CP[as.vector(mycl.C2CP)]

# assign your favorite heatmap color scheme. Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]). Type demo.col(20) to see more color schemes.
myheatcol <- colorRampPalette(colors=c("yellow","white","blue"))(100)
# plot the hclust results as a heatmap
heatmap.2(diffSets.C2CP,
          Rowv=as.dendrogram(hr.C2CP),
          Colv=NA,
          col=myheatcol, scale="row",
          density.info="none", trace="none",
          cexRow=0.9, cexCol=1, margins=c(10,14))
