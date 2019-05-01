# read in Kallisto quants ----
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(biomaRt) # provides access to a wealth of annotation info
targets <- read_tsv("studyDesign.txt")# read in your study design
path <- file.path(targets$sample, "abundance.h5") # set file paths to your mapped data
targets <- mutate(targets, path) # add paths to your study design (only necessary for Sleuth)
Hs.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl") # select 'mart' from biomaRt for annotations
Tx <- getBM(attributes=c('ensembl_transcript_id_version', # get gene symbols for each transcript ID
                         'external_gene_name'),
            mart = Hs.anno)
Tx <- as_tibble(Tx) # convert this annotation mapping file to a tibble (the tidyverse version of a dataframe)
Txi_gene <- tximport(path, #reading kallisto data into R
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, 
                     countsFromAbundance = "lengthScaledTPM")
myCPM <- as_tibble(Txi_gene$abundance, rownames = "geneSymbol") # these are you counts after adjusting for transcript length
myCounts <- as_tibble(Txi_gene$counts, rownames = "geneSymbol") # these are your transcript per million (TPM) values, or counts per million (CPM) if you collapsed data to gene level

# Filter and normalize ----
library(RColorBrewer) 
library(reshape2) 
library(genefilter)
library(edgeR) 
library(matrixStats)
library(hrbrthemes)
groups <- targets$condition
groups <- factor(groups)
sampleLabels <- targets$sample
myDGEList <- DGEList(Txi_gene$counts)
save(myDGEList, file = "myDGEList")
load(file = "myDGEList")
log2.cpm <- cpm(myDGEList, log=TRUE)
nsamples <- ncol(log2.cpm)
myColors <- brewer.pal(nsamples, "Paired")
log2.cpm.df <- as_tibble(log2.cpm)
colnames(log2.cpm.df) <- sampleLabels
log2.cpm.df.melt <- melt(log2.cpm.df)

ggplot(log2.cpm.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.y = "median", geom = "point", shape = 124, size = 6, color = "black", show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  coord_flip() +
  theme_ipsum_rc() 

table(rowSums(myDGEList$counts==0)==10)
cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=5 #user defined
myDGEList.filtered <- myDGEList[keepers,]
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm)
colnames(log2.cpm.filtered.norm.df) <- sampleLabels
log2.cpm.filtered.norm.df.melt <- melt(log2.cpm.filtered.norm.df)

ggplot(log2.cpm.filtered.norm.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.y = "median", geom = "point", shape = 124, size = 6, color = "black", show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  coord_flip() +
  theme_ipsum_rc() 

# plot PCA ----
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
x <- pca.res$rotation 
pc.var<-pca.res$sdev^2
pc.per<-round(pc.var/sum(pc.var)*100, 1)

pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df, aes(x=PC1, y=PC2, color=groups)) +
  geom_point(size=4) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_ipsum_rc()

# differential expression ----
library(limma) #powerful package for differential gene expression using linear modeling

design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)
fit <- lmFit(v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts(infection = infected - control,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits <- as_tibble(myTopHits, rownames = "geneSymbol")
ggplot(myTopHits, aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneSymbol))) +
           geom_point(size=2) +
           ylim(-0.5,12) +
           geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
           geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
           geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
           labs(title="Volcano plot",
                subtitle = "impact of L. braziliensis infection on human lncRNAs",
                caption=paste0("produced on ", Sys.time())) +
           theme_ipsum_rc()

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=2)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneSymbol")

# module identification ----
library(gt)
library(gplots) #the heatmap2 function in this package is a primary tool for making heatmaps

myheatcolors2 <- colorRampPalette(colors=c("yellow","white","blue"))(100)
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") #cluster rows by pearson correlation
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=NA,
          RowSideColors=module.color,
          col=myheatcolors2, scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20)) 

# view modules of co-regulated genes
# view your color assignments for the different clusters
names(module.color) <- names(module.assign) 
barplot(rep(10, max(module.assign)),
        col=unique(module.color[clustRows$labels[clustRows$order]]), 
        horiz=T, names=unique(module.assign[clustRows$order]))

#choose a cluster(s) of interest by selecting the corresponding number based on the previous graph
module.up <- 1 
module.down <- 2
myModule.up <- diffGenes[names(module.assign[module.assign%in%module.up]),] 
myModule.down <- diffGenes[names(module.assign[module.assign%in%module.down]),] 
hrsub.up <- hclust(as.dist(1-cor(t(myModule.up), method="pearson")), method="complete") 
hrsub.down <- hclust(as.dist(1-cor(t(myModule.down), method="pearson")), method="complete") 

#prints out genes in the order you see them in the cluster
moduleSymbols.up <- data.frame(Labels=rev(hrsub.up$labels[hrsub.up$order]))
write_tsv(moduleSymbols.up, "moduleSymbols.up.txt")
moduleSymbols.up <- as.vector(t(moduleSymbols.up))

moduleSymbols.down <- data.frame(Labels=rev(hrsub.down$labels[hrsub.down$order]))
moduleSymbols.down <- as.vector(t(moduleSymbols.down))

moduleData.up <- diffGenes[moduleSymbols.up,]
moduleData.up.df <- as_tibble(moduleData.up, rownames = "geneSymbol")

moduleData.down <- diffGenes[moduleSymbols.down,]
moduleData.down.df <- as_tibble(moduleData.down, rownames = "geneSymbol")


# Use dplyr 'mutate' function to add new columns based on existing data
moduleData.up.df <- mutate(moduleData.up.df,
                    control_avg = (control_04 + control_05 + control_06 + control_07 + control_08)/5, 
                    infected_avg = (patient_26912 + patient_26968 + patient_26969 + patient_26971 + patient_26972)/5,
                    LogFC = (infected_avg - control_avg)) %>%
  mutate_if(is.numeric, round, 2)

#now look at this modified data table
moduleData.up.df

# Make a table with the top 5 ncRNAs induced by infection
mydata.up.sort <- moduleData.up.df %>%
  dplyr::arrange(desc(LogFC)) %>% #note that this is the first time you've seen the 'pipe' operator
  dplyr::select(geneSymbol, control_avg, infected_avg, LogFC) %>%
  dplyr::top_n(5)

gt(mydata.up.sort %>%
     tab_header(
       title = md("**Top 5 ncRNAs induced by L. braziensis**")
     ))

# # Repeat the above, but with the down regulated module
moduleData.down.df <- mutate(moduleData.down.df,
                           control_avg = (control_04 + control_05 + control_06 + control_07 + control_08)/5, 
                           infected_avg = (patient_26912 + patient_26968 + patient_26969 + patient_26971 + patient_26972)/5,
                           LogFC = (infected_avg - control_avg)) %>%
  mutate_if(is.numeric, round, 2)

#now look at this modified data table
moduleData.down.df

mydata.down.sort <- moduleData.down.df %>%
  dplyr::arrange(LogFC) %>% #note that this is the first time you've seen the 'pipe' operator
  dplyr::select(geneSymbol, control_avg, infected_avg, LogFC) %>%
  dplyr::top_n(5)

gt(mydata.down.sort %>%
     tab_header(
       title = md("**Top 5 ncRNAs repressed by L. braziensis**")
     ))

