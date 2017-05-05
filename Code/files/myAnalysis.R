#this is where you'll build your analysis script throughout the workshop
##############################################################################################################################
#Step1_preprocessRNAseq
#This script carries out the steps involved in analysis of RNAseq data.  
#depending on your data and interests, different parts of this script may not apply to you
##############################################################################################################################
#begin by loading the packages required for RNAseq data
library(Rsubread)
library(limma)
library(edgeR)
library(ShortRead)
options(digits=2)
#library(Biostrings)

#read in your study design file
targets <- readTargets("Igor_Ripk3_BMMC_studyDesign_RNAseq.txt", row.names=NULL)
targets
targets.mod <- targets[1:12,]
groups <- factor(paste(targets.mod$genotype, targets.mod$treatment, sep="."))
#create some more human-readable labels for your samples using the info in this file
sampleLabels <- paste(targets.mod$genotype, targets.mod$treatment, targets.mod$rep, sep=".")

#set-up your experimental design
design <- model.matrix(~0+groups)
colnames(design) <- levels(groups)
design

# ##############################################################################################################################
# #you can check read quality using shortRead package
# #but I usually find it is better to do this on the sequencer or Illumina's BaseSpace website
# ##############################################################################################################################
# myFastq <- targets$fastq
# #collecting statistics over the files
# qaSummary <- qa("B6-WT-untreat-rep2_S2_mergedLanes_read1.fastq", type="fastq")
# #create and view a report
# browseURL(report(qaSummary))

##############################################################################################################################
#build index from your reference genome (expect this to take about 20 min on 8G RAM for mouse genome)
#you must have already downloaded the fasta file for your genome of interest and have it in the working directory
#this only needs to be done once, then index can be reused for future alignments
##############################################################################################################################
buildindex(basename="mouse",reference="Mus_musculus.GRCm38.dna.primary_assembly.fa")

##############################################################################################################################
#align your reads (in the fastq files) to your indexed reference genome that you created above
#expect this to take about 45min for a single fastq file containing 25 million reads
#the output from this is a .bam file for each of your original fastq files
##############################################################################################################################
reads1 <- targets$fastq[9]
reads2 <- targets$fastq[21] 
align(index="mouse", readfile1=reads1, readfile2=reads2, input_format="gzFASTQ",output_format="BAM",
      output_file="alignmentResultsPE_sample9.BAM", tieBreakHamming=TRUE,unique=TRUE,indels=5, nthreads=8)

##############################################################################################################################
#use the 'featureCounts' function to summarize read counts to genomic features (exons, genes, etc)
#will take about 1-2min per .bam file.
#for total transcriptome data summarized to mouse/human .gtf, expect about 50-60% of reads summarize to genes (rest is non-coding)
##############################################################################################################################
#read in text file with bam file names
MyBAM <- read.delim("BAMfiles_names.txt", header=T)
MyBAM <- as.character(MyBAM[,1])
#summarize aligned reads to genomic features (i.e. exons)
fc <- featureCounts(files=MyBAM, annot.ext="Mus_musculus.GRCm38.79.gtf", isGTFAnnotationFile=TRUE, GTF.featureType = "exon",
                    GTF.attrType="gene_id", useMetaFeatures=TRUE, isPairedEnd=TRUE, requireBothEndsMapped=TRUE, strandSpecific=2, nthreads=8)
#use the 'DGEList' function from EdgeR to make a 'digital gene expression list' object
DGEList <- DGEList(counts=fc$counts, genes=fc$annotation)
load("DGElist")
#retrieve all your gene/transcript identifiers from this DGEList object
myEnsemblIDs <- DGEList$genes$GeneID

##############################################################################################################################
#Normalize unfiltered data using 'voom' function in Limma package
#This will normalize based on the mean-variance relationship
#will also generate the log2 of counts per million based on the size of each library (also a form of normalization)
##############################################################################################################################
normData.unfiltered <- voom(DGEList, design, plot=TRUE)
exprs.unfiltered <- normData.unfiltered$E
exprs.matrix.unfiltered <- as.matrix(exprs.unfiltered)
#note that because you're now working with Log2 CPM, much of your data will be negative number (log2 of number smaller than 1 is negative) 
head(exprs.matrix.unfiltered)

#if you need RPKM for your unfiltered, they can generated as follows
#Although RPKM are commonly used, not really necessary since you don't care to compare two different genes within a sample
rpkm.unfiltered <- rpkm(DGEList, DGEList$genes$Length)
tail(rpkm)

##############################################################################################################################
#Filtering your dataset
#Only keep in the analysis those genes which had >10 reads per million mapped reads in at least two libraries.
#Need to remove duplicate gene symbols
##############################################################################################################################
cpm.matrix.filtered <- rowSums(cpm(DGEList) > 10) >= 2
DGEList.filtered <- DGEList[cpm.matrix.filtered,]
dim(DGEList.filtered)
rpkm.filtered <- rpkm(DGEList.filtered, DGEList.filtered$genes$Length) #if you prefer, can use 'cpm' instead of 'rpkm' here
# "+1" because log2 of negative and 0 will give you infinity and you want to avoid that
rpkm.filtered <- log2(rpkm.filtered + 1)

#Use Voom again to normalize this filtered dataset
normData <- voom(DGEList.filtered, design, plot=TRUE)
exprs <- normData$E
exprs.matrix.filtered <- as.matrix(exprs)
head(exprs.matrix.filtered)

##############################################################################################################################
#annotate your normalized data using the organism-specific database package
##############################################################################################################################
library(org.Mm.eg.db)
library(AnnotationDbi)
#ls("package:org.Mm.eg.db")
#If we want to know what kinds of data are retriveable via the 'select' command, look at the columns of the annotation database
columns(org.Mm.eg.db)
#If we want to know what kinds of fields we could potentially use as keys to query the database, use the 'keytypes' command
keytypes(org.Mm.eg.db)
#transform your identifiers to entrezIDs
myAnnot.unfiltered <- AnnotationDbi::select(org.Mm.eg.db, keys=rownames(exprs.matrix.unfiltered), keytype="ENSEMBL", columns=c("ENTREZID", "REFSEQ", "UNIGENE", "GENENAME", "SYMBOL"))
myAnnot.filtered <- AnnotationDbi::select(org.Mm.eg.db, keys=rownames(exprs.matrix.filtered), keytype="ENSEMBL", columns=c("ENTREZID", "REFSEQ", "UNIGENE", "GENENAME", "SYMBOL"))
resultTable.unfiltered <- merge(myAnnot.unfiltered, exprs.matrix.unfiltered, by.x="ENSEMBL", by.y=0)
resultTable.filtered <- merge(myAnnot.filtered, exprs.matrix.filtered, by.x="ENSEMBL", by.y=0)
head(resultTable.unfiltered)
#add more appropriate sample names as column headers
colnames(resultTable.unfiltered) <- c("Ensembl", "Entrez", "Refseq", "Unigene", "Name", "Symbol", sampleLabels)
colnames(resultTable.filtered) <- c("Ensembl", "Entrez", "Refseq", "Unigene", "Name", "Symbol", sampleLabels)
#now write these annotated datasets out
write.table(resultTable.unfiltered, "normalizedUnfilteredRNAseq.txt", sep="\t", quote=FALSE)
write.table(resultTable.filtered, "normalizedFilteredRNAseq.txt", sep="\t", quote=FALSE)
head(resultTable.unfiltered)
head(resultTable.filtered)

#end of Step1_preprocessRNAseq

###########################################################
#Start of Collapse Function
#I prefer to use a table of data with a non-redundant set of gene identifiers as input for this script
#this keeps my heatmaps and lists of DEGS as simple of as possible
#if you are coming to this script with array data, you've already reduced your list in this way
#RNAseq data is a bit more complicated (many transcripts, but the same gene), 
#but you achieve the same level of reduction by filtering based on annotation data
#################################################################

#additional filtering based on annotation
#if you have array data, you can skip this section
#################################################################

#our initial filtering of RNAseq data was based on solely cpm, which only removes relatively lowly expressed genes
#now that you have annotation info, you can further filter to reduce to single line of data per gene symbo using one of two options
#first, check out how many unqiue genes are represented in your data
dupFiltered <- unique(resultTable.filtered$Symbol)
#use collapseRows function from WGCNA package to collapse your dataset
library(WGCNA)
#pull your rownames and unique identifiers
myIDs <- rownames(resultTable.filtered)
#retrieve your gene symbols from your data
mySymbols <- resultTable.filtered[,6]
#remove all annotation columns so you're left with only numeric data
resultTable <- resultTable.filtered[,-1:-6]
myCollapsed <- collapseRows(resultTable, mySymbols, myIDs, method = "MaxMean")
myCollapsed <- myCollapsed$datETcollapsed
write.table(myCollapsed, "Collapsed.txt", sep="\t", quote=FALSE)
#now that the matrix is collapsed to give non-redudant list of genes,
#you could set the symbols to be the row names and move on

#for GSEA tab delimited file. Clean data (remove unnecessary columns and make lower case upper case)
myCollapsedRNAseq <- read.delim("Collapsed.txt", header=T, sep="\t")
head(myCollapsedRNAseq)
#need to convert gene symbols to uppercase
mySymbolsforGSEA <- toupper(row.names(myCollapsedRNAseq))
head(mySymbolsforGSEA)
#need to add column back to data
myCollapsedRNAseq <- cbind(mySymbolsforGSEA, myCollapsedRNAseq)
write.table(myCollapsedRNAseq, "myCollapsedRNAseq.txt", sep="\t", quote=FALSE)

#to make grp file from other published arrays for GSEA and make lower case upper case)
FontanaJunB <- read.delim("FontanaJunB.txt", header=T, sep="\t")
head(FontanaJunB)
#need to convert gene symbols to uppercase
FontanaJunB <- toupper(FontanaJunB[,1])
dim(FontanaJunB)
head(FontanaJunB)
write.table(FontanaJunB, "FontanaJunB.txt", sep="\t", quote=FALSE)
##############################################################################################

#Step2_dataExploration
#goal of this script is to using multivariate statisical approaches to explore the structure of your data
#begin by taking a look at the text expression matrix you created at the end of the last class
#head(exprs.matrix.filtered) # don't need this because I have myCollapsed data

###############################################################################################
# set up your experimental design by reading in a targets file that explains treatments, conditions, etc
###############################################################################################
library(limma)
#read in a tab-delimited "targets" file with the study design
targets <- readTargets("Igor_Ripk3_BMMC_studyDesign_RNAseq.txt", sep="\t")
targets
myGroups <- factor(paste(targets$genotype, targets$treatment, sep="."))
myGroups
design <- model.matrix(~0+myGroups)
colnames(design) <- levels(myGroups)
design


###############################################################################################
#carry out hierarchical clustering on filtered data
###############################################################################################
#make some sample labels 
sampleLabels <- paste(targets$genotype, targets$treatment, sep=".")
distance <- dist(t(exprs.matrix.filtered),method="maximum")
clusters <- hclust(distance, method = "complete") 
#plot cluster dendrogram
plot(clusters, label = sampleLabels, hang = -1)


###############################################################################################
#Principal component analysis of the filtered data matrix T=transpose
###############################################################################################
pca.res <- prcomp(t(exprs.matrix.filtered), scale.=F, retx=T)
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
head(pca.res$rotation) #$rotation shows you how much each GENE influenced each PC (callled 'eigenvalues', or loadings)
head(pca.res$x) #$x shows you how much each SAMPLE influenced each PC (called 'scores')
plot(pca.res, las=1) #las sets labels horizontal
pc.var<-pca.res$sdev^2 #sdev^2 gives you the eigenvalues
pc.per<-round(pc.var/sum(pc.var)*100, 1) #percentage variance
pc.per

#make some graphs to visualize your PCA result
library(ggplot2)
#lets first plot any two PCs against each other
#turn your scores for each gene into a data frame
data.frame <- as.data.frame(pca.res$x)
ggplot(data.frame, aes(x=PC1, y=PC2, colour=factor(myGroups))) +
  geom_point(size=5) +
  theme(legend.position="right")

#create a 'small multiples' chart to look at impact of each variable on each pricipal component
library(reshape2)
melted <- cbind(myGroups, melt(pca.res$x[,1:3]))
#look at your 'melted' data
ggplot(melted) +
  geom_bar(aes(x=Var1, y=value, fill=myGroups), stat="identity") +
  facet_wrap(~Var2)
#if you have a batch effect can use combat function to take care of that

###################################################################################################
#########################################

### Step3_dataWrangling

#this script walks thorough some basic data wrangling for organizing expression data spreadsheets and ends with
#how to create publication quality graphics from transcriptomic data generated (regardless of platform used)
#to start this script you need a file with all your expression data and some non-redundant identifiers as row names (usually gene symbols)
#you also need a study design file

#for creating graphics, we'll use the ggplot2 and ggvis packages which employ a 'grammar of graphics' approach
#load the packages
library(ggplot2)
library(dplyr)
library(ggvis)

#read in your data from a text file that contains genes symbols as rows, and samples as columns. 
myCollapsedData <- read.delim("Collapsed.txt", header=TRUE)
head(myCollapsedData)
#subset data to get rid of column "2" = -2
#myData <- myData[,-2]
#make column 1 the row names of myData
#row.names(myData) <- myData[,1]
#myData <- myData[,-1]
#head(myData)

# [This is probably not necessary since myCollapsedData already has column headers, but for array data use this] column headers are a bit cumbersome, so we'll change these to something more human-readable
#targets <- read.delim("Igor_Ripk3_BMMC_studyDesign_RNAseq.txt", sep="\t", stringsAsFactors = FALSE)
#sampleLabels <- as.character(paste(targets$cellType, targets$treatment, targets$phenotype, targets$rep, sep="."))
#colnames(myData) <- sampleLabels
#head(myData)
geneSymbols <- row.names(myCollapsedData)
#myCollapsedData.df <- as.data.frame(myCollapsedData)
#head(myCollapsedData.df)

#use the dplyr 'mutate' command to get averages and fold changes for all your replicates, dplyr only works on data frames and row names are not considered part of the table
myCollapsedData <- mutate(myCollapsedData,
                      B6.untreated.AVG = (B6.untreated.1 + B6.untreated.2)/2,
                      Ripk3.untreated.AVG = (Ripk3.untreated.1 + Ripk3.untreated.2)/2,
                      Ripk3_Casp8.untreated.AVG = (Ripk3_Casp8.untreated.1 + Ripk3_Casp8.untreated.2)/2,
                      B6.LPS_6hr.AVG = (B6.LPS_6hr.1 + B6.LPS_6hr.2)/2,
                      Ripk3.LPS_6hr.AVG = (Ripk3.LPS_6hr.1 + Ripk3.LPS_6hr.2)/2,
                      Ripk3_Casp8.LPS_6hr.AVG = (Ripk3_Casp8.LPS_6hr.1 + Ripk3_Casp8.LPS_6hr.2)/2,
                      LogFC.Ripk3.LPS_6hr.vs.B6.LPS_6hr = (Ripk3.LPS_6hr.AVG - B6.LPS_6hr.AVG),
                      LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr = (Ripk3_Casp8.LPS_6hr.AVG - Ripk3.LPS_6hr.AVG),
                      LogFC.Ripk3_Casp8.LPS_6hr.vs.B6.LPS_6hr = (Ripk3_Casp8.LPS_6hr.AVG - B6.LPS_6hr.AVG),
                      LogFC.Ripk3_Casp8.untreated.vs.Ripk3.untreated = (Ripk3_Casp8.untreated.AVG - Ripk3.untreated.AVG),
                      LogFC.Ripk3.untreated.vs.B6.untreated = (Ripk3.untreated.AVG - B6.untreated.AVG),
                      LogFC.Ripk3_Casp8.untreated.vs.B6.untreated = (Ripk3_Casp8.untreated.AVG - B6.untreated.AVG),
                      geneSymbols=geneSymbols) # why do we have a column at the end with gene Symbols??

#now look at this modified data table
row.names(myCollapsedData) <- geneSymbols
head(myCollapsedData)
dim(myCollapsedData)
#capture the symbols from the rownames and set to an object called 'symbols'
symbols <- rownames(myCollapsedData)
#now get rid of rownames
rownames(myCollapsedData) <- NULL
#now edit your data to contain a first column with symbols
myCollapsedData <- cbind(symbols, myCollapsedData)
head(myCollapsedData)
#now save this new data object to replace the old one
save(myCollapsedData, file="myCollapsedData")


#use dplyr "arrange" and "select" functions to sort by LogFC column of interest (arrange)
#and then display only the columns of interest (select) to see the most differentially expressed genes
myCollapsedData.sort <- myCollapsedData %>%
  arrange(desc(LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr)) %>%
  select(geneSymbols, LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr) 
head(myCollapsedData.sort)

#use dplyr "filter" and "select" functions to pick out genes of interest (filter)
#ways to tweek the 'select' function
#use : between two column names to select all columns between
#use 'contains', 'starts_with' or 'ends_with' to modify how you select
#can refer to columns using exact name or numerical indicator
myCollapsedData.filter <- myCollapsedData %>%
  filter(geneSymbols=="Irf1" | geneSymbols=="Junb") %>%
  select(geneSymbols, LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr) 
head(myCollapsedData.filter)

myCollapsedData.filter <- myCollapsedData %>%
  filter(geneSymbols=="Irf1" | geneSymbols=="Junb") %>%
  select(geneSymbols, LogFC.Ripk3_Casp8.LPS_6hr.vs.B6.LPS_6hr, LogFC.Ripk3.LPS_6hr.vs.B6.LPS_6hr, LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr ) 
head(myCollapsedData.filter)

#another example using grep to match patterns
myCollapsedData.filter <- myCollapsedData %>%
  #filter(grepl('Ccl', geneSymbols)) %>%
  filter(grepl('Il', geneSymbols)) %>%
  select(geneSymbols, LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr) 
head(myCollapsedData.filter)

#first reorder and clean up the data you want to graph
#row.names(myCollapsedData.filter) <- myCollapsedData.filter[,1]
#myCollapsedData.filter <- select(myCollapsedData.filter, -geneSymbols)
#myCollapsedData.filter.transpose <- as.data.frame(t(myCollapsedData.filter))
#cellType <- c("B1a", "B1b", "Mac")
#myCollapsedData.filter.transpose <- mutate(myData.filter.transpose, cellType, 
                                  phenotype=c("naive", "naive","naive"))
#myCollapsedData.filter.transpose

#another way to set up for graphing...don't filter before transposing
#myCollapsedData.graph <- select(myCollapsedData,-geneSymbols)
#head(myCollapsedData.graph)
#myCollapsedData.graph.transpose <- as.data.frame(t(myCollapsedData.graph))
#dim(myCollapsedData.graph.transpose)

#plot a simple bar graph - y=column of interest, x= row of interest from filtered table
ggplot(myCollapsedData.filter, aes(y=LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr)) +
  geom_bar(aes(x=geneSymbols, fill=LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr), stat="identity") +
  theme(axis.text.x=element_text(angle=-45))

#create a basic scatterplot using ggplot
ggplot(myCollapsedData, aes(x=Ripk3.LPS_6hr.AVG, y=Ripk3_Casp8.LPS_6hr.AVG)) +
  geom_point(shape=1) +
  geom_point(size=4)

# ##Volcano plots
# ##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
# gene_list$threshold = as.factor(abs(gene_list$logFC) > 2 & gene_list$P.Value < 0.05/no_of_genes)
# 
# ##Construct the volcano plot
# ggplot(data=myData, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
#   geom_point(alpha=0.4, size=1.75) +
#   opts(legend.position = "none") +
#   xlim(c(-10, 10)) + ylim(c(0, 15)) +
#   xlab("log2 fold change") + ylab("-log10 p-value")

#define a tooltip that shows gene symbol and Log2 expression data when you mouse over each data point in the plot
tooltip <- function(data, ...) {
  paste0("<b>","Symbol: ", data$geneSymbols, "</b><br>",
         "Ripk3.LPS_6hr.AVG: ", data$Ripk3.LPS_6hr.AVG, "<br>",
         "Ripk3_Casp8.LPS_6hr.AVG: ", data$Ripk3_Casp8.LPS_6hr.AVG)
}

#plot the interactive graphic
myCollapsedData %>% 
  ggvis(x= ~Ripk3.LPS_6hr.AVG, y= ~Ripk3_Casp8.LPS_6hr.AVG, key := ~geneSymbols) %>% 
  layer_points(fill = ~LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr) %>%
  add_tooltip(tooltip)

################################################################

#Step 4

#first part of this script is in Step 1 
#the goal of this script is to identify differentially expressed genes (DEG) 
#you should already know what pairwise comparisons are most important to you

#I prefer to use a table of data with a non-redundant set of gene identifiers as input for this script
#this keeps my heatmaps and lists of DEGS as simple of as possible
#if you are coming to this script with array data, you've already reduced your list in this way
#RNAseq data is a bit more complicated (many transcripts, but the same gene), 
#but you achieve the same level of reduction by filtering based on annotation data
#################################################################
#additional filtering based on annotation
#if you have array data, you can skip this section
#################################################################
#our initial filtering of RNAseq data was based on solely cpm, which only removes relatively lowly expressed genes
#now that you have annotation info, you can further filter to reduce to single line of data per gene symbol
#first, check out how many unqiue genes are represented in your data
dupFiltered <- unique(resultTable.filtered$Symbol)
#use collapseRows function from WGCNA package to collapse your dataset
library(WGCNA)
#pull your rownames and unique identifiers
myIDs <- rownames(resultTable.filtered)
#retrieve your gene symbols from your data
mySymbols <- resultTable.filtered[,6]
#remove all annotation columns so you're left with only numeric data because you need to calculate the mean
resultTable <- resultTable.filtered[,-1:-6]
# function(data table, the group that is redundant, unique identifier, method of collapsing)
myCollapsed <- collapseRows(resultTable, mySymbols, myIDs, method = "MaxMean")
myCollapsed <- myCollapsed$datETcollapsed
write.table(myCollapsed, "Collapsed.txt", sep="\t", quote=FALSE)
#now that the matrix is collapsed to give non-redudant list of genes,
#you could set the symbols to be the row names and move on

#################################################################################################################
#if you have no biological replicates, you will not be able to leverage statistical tools to identify DE genes
#Instead, you will ONLY rely on fold changes
#if you DO have replicates, skip this section and proceed to the next part 
#################################################################################################################
#use the dplyr 'filter' command to capture all the genes that are up/down regulated x-fold in n conditions
#in this case, 'myData' is a dataframe that you generated with Log2 expression and annotation
#myData.filter <- myData %>%
  #filter((abs(Ecdysone.vs.PBS_18hr_gut) >= 1) | (abs(Ecdysone.vs.PBS_5hr_gut) >= 1)) %>%
  #select(geneID, Ecdysone.vs.PBS_5hr_carcass, Ecdysone.vs.PBS_18hr_carcass, Ecdysone.vs.PBS_5hr_gut, Ecdysone.vs.PBS_18hr_gut)
#head(myData.filter)


###############################################################################################
# use Limma to find differentially expressed genes between two or more conditions
###############################################################################################
# fit the linear model to your filtered expression data in matrix format only samples and rownames no AVG or FC
library(limma)
rownames(myCollapsed)
dim(myCollapsed)
fit <- lmFit(myCollapsed, design)
#feed in data from the differential gene expression of UT vs LPS
rownames(diffDataUT.LPS)
fitUTvLPS <- lmFit(diffDataUT.LPS, design)
#add annotation into your linear model fit
#don't really need to do this if you have RNAseq data
#library(annotate)
#fit$genes$Symbol <- getSYMBOL(probeList, "lumiMouseAll.db")
#fit$genes$Entrez <- getEG(probeList, "lumiMouseAll.db")

# set up a contrast matrix based on the pairwise comparisons of interest
contrast.matrix.LPS <- makeContrasts(R3C8.v.R3 = Ripk3_Casp8.LPS_6hr - Ripk3.LPS_6hr, 
                                     R3.v.B6 = Ripk3.LPS_6hr - B6.LPS_6hr, 
                                     R3C8.v.B6 = Ripk3_Casp8.LPS_6hr - B6.LPS_6hr, levels=design)
contrast.matrix.untreated <- makeContrasts(R3C8.v.R3.untreated = Ripk3_Casp8.untreated - Ripk3.untreated, 
                                           R3.v.B6.untreated = Ripk3.untreated - B6.untreated, R3C8.v.B6.untreated = Ripk3_Casp8.untreated - B6.untreated, levels=design)

### set up a contrast matrix comparing UT and LPS
contrast.matrix.UT.LPS <- makeContrasts(B6UT.v.B6LPS = B6.LPS_6hr - B6.untreated, R3UT.v.R3LPS = Ripk3.LPS_6hr - Ripk3.untreated, R3C8UT.v.R3C8LPS = Ripk3_Casp8.LPS_6hr - Ripk3_Casp8.untreated, levels=design)

# check each contrast matrix
contrast.matrix.LPS
contrast.matrix.untreated
contrast.matrix.UT.LPS

# extract the linear model fit for the contrast matrix that you just defined above
fitsLPS <- contrasts.fit(fit, contrast.matrix.LPS)
fitsUntreated <- contrasts.fit(fit, contrast.matrix.untreated)
fitsUT.LPS <- contrasts.fit(fit, contrast.matrix.UT.LPS)

# Now do the same for the UT vs LPS Diff Expr Genes
fits.of.LPS <- contrasts.fit(fitUTvLPS, contrast.matrix.LPS)

#get bayesian stats for your linear model fit
ebFitLPS <- eBayes(fitsLPS)
ebFitUntreated <- eBayes(fitsUntreated)
ebFitUT.LPS <- eBayes(fitsUT.LPS)

#Now get the bayesian stats for the UTvsLPS data set
ebFit.of.LPS <- eBayes(fits.of.LPS)
###############################################################################################
# use the topTable and decideTests functions to see the differentially expressed genes
###############################################################################################

# use topTable function to take a look at the hits (BH=benjamini hochberg?)
#coef = which of the contrast matrix columns are you interested in?
#number = how many genes are you interested in?
myTopHitsLPS <- topTable(ebFitLPS, adjust ="BH", coef=3, number=20, sort.by="logFC")
myTopHitsLPS
myTopHitsUntreated <- topTable(ebFitUntreated, adjust ="BH", coef=3, number=20, sort.by="logFC")
myTopHitsUntreated
myTopHitsUT.LPS <- topTable(ebFitUT.LPS, adjust ="BH", coef=2, number=20, sort.by="logFC")
myTopHitsUT.LPS

myTopHits.of.LPS <- topTable(ebFit.of.LPS, adjust ="BH", coef=3, number=20, sort.by="logFC")
myTopHits.of.LPS

# use the 'decideTests' function to show Venn diagram for all diffexp genes for up to three comparisons
# 0 = not differentially expressed, 1 = up, -1 = down. Now show me all the 1s. lfc = log fold change
resultsLPS <- decideTests(ebFitLPS, method="global", adjust.method="BH", p.value=0.05, lfc=0.59)
#stats <- write.fit(ebFitLPS)
vennDiagram(resultsLPS, include="both") #all pairwise comparisons on a B6 background
resultsUntreated <- decideTests(ebFitUntreated, method="global", adjust.method="BH", p.value=0.05, lfc=0.59)
vennDiagram(resultsUntreated, include="both", cex=1) #all pairwise comparisons on a B6 background

resultsUT.LPS <- decideTests(ebFitUT.LPS, method="global", adjust.method="BH", p.value=0.05, lfc=0.59)
vennDiagram(resultsUT.LPS, include="both", cex = 1)

# take a look at what the results of decideTests looks like
resultsLPS
resultsUT.LPS

# now pull out probeIDs from selected regions of the Venn diagram.  In this case, I want all genes in the venn. The columns refer to the contrast matrix
diffProbesLPS <- which(resultsLPS[,1] !=0 | resultsLPS[,2] !=0 | resultsLPS[,3] !=0)
class(diffProbesLPS)
diffProbesUntreated <- which(resultsUntreated[,1] !=0 | resultsUntreated[,2] !=0 | resultsUntreated[,3] !=0)
diffProbesUT.LPS <- which(resultsUT.LPS[,1] !=0 | resultsUT.LPS[,2] !=0 | resultsUT.LPS[,3] !=0)

#before pulling out expression data for differentially expressed genes, convert matrix to eset with annotation
library(Biobase)
myEset.ALL <- new("ExpressionSet", exprs = myCollapsed)

# retrieve expression data for the probes from above
diffData.LPS <- myEset.ALL[resultsLPS[,1] !=0 | resultsLPS[,2] !=0 | resultsLPS[,3] !=0]
diffData.Untreated <- myEset.ALL[resultsUntreated[,1] !=0 | resultsUntreated[,2] !=0 | resultsUntreated[,3] !=0]
diffDataUT.LPS <- myEset.ALL[resultsUT.LPS[,1] !=0 | resultsUT.LPS[,2] !=0 | resultsUT.LPS[,3] !=0]

#pull the expression data back out of the eset object
diffData.LPS <- exprs(diffData.LPS)
diffData.Untreated <- exprs(diffData.Untreated)
diffDataUT.LPS <- exprs(diffDataUT.LPS)

#combine probeIDs, gene symbols and expression data for differentially expressed genes into one file
write.table(cbind(diffProbesLPS, diffData.LPS),"DiffGenesLPS.xls", sep="\t", quote=FALSE)
write.table(cbind(diffProbesUntreated, diffData.Untreated),"DiffGenesUntreated.xls", sep="\t", quote=FALSE)
## why do we include diffProbes as one of the columns here? cause that's the one with all the weird numbers, no?

# take a look at each expression matrix
dim(diffData.LPS)
dim(diffData.Untreated)
dim(diffDataUT.LPS)

#############################################################################################################

### Step3b_dataWrangling_for Diff expressed genes

#this script walks thorough some basic data wrangling for organizing expression data spreadsheets and ends with
#how to create publication quality graphics from transcriptomic data generated (regardless of platform used)
#to start this script you need a file with all your expression data and some non-redundant identifiers as row names (usually gene symbols)
#you also need a study design file

#for creating graphics, we'll use the ggplot2 and ggvis packages which employ a 'grammar of graphics' approach
#load the packages
library(ggplot2)
library(dplyr)
library(ggvis)

#read in your data from a text file that contains genes symbols as rows, and samples as columns. 
myDiffGenes.LPS <- read.delim("LPS.Cluster3.1.txt", header=TRUE)
head(myDiffGenes.LPS)
#myDiffGenes.LPS <- myDiffGenes.LPS[,-1]
#subset data to get rid of column "2" = -2
#myData <- myData[,-2]
#make column 1 the row names of myData
#row.names(myData) <- myData[,1]
#myData <- myData[,-1]
#head(myData)

# [This is probably not necessary since myCollapsedData already has column headers, but for array data use this] column headers are a bit cumbersome, so we'll change these to something more human-readable
#targets <- read.delim("Igor_Ripk3_BMMC_studyDesign_RNAseq.txt", sep="\t", stringsAsFactors = FALSE)
#sampleLabels <- as.character(paste(targets$cellType, targets$treatment, targets$phenotype, targets$rep, sep="."))
#colnames(myData) <- sampleLabels
#head(myData)
geneSymbols <- (myDiffGenes.LPS[,1])
#myCollapsedData.df <- as.data.frame(myCollapsedData)
#head(myCollapsedData.df)

#use the dplyr 'mutate' command to get averages and fold changes for all your replicates, dplyr only works on data frames and row names are not considered part of the table
myDiffGenes.LPS <- mutate(myDiffGenes.LPS,
                          B6.untreated.AVG = (B6.untreated.1 + B6.untreated.2)/2,
                          Ripk3.untreated.AVG = (Ripk3.untreated.1 + Ripk3.untreated.2)/2,
                          Ripk3_Casp8.untreated.AVG = (Ripk3_Casp8.untreated.1 + Ripk3_Casp8.untreated.2)/2,
                          B6.LPS_6hr.AVG = (B6.LPS_6hr.1 + B6.LPS_6hr.2)/2,
                          Ripk3.LPS_6hr.AVG = (Ripk3.LPS_6hr.1 + Ripk3.LPS_6hr.2)/2,
                          Ripk3_Casp8.LPS_6hr.AVG = (Ripk3_Casp8.LPS_6hr.1 + Ripk3_Casp8.LPS_6hr.2)/2,
                          LogFC.Ripk3.LPS_6hr.vs.B6.LPS_6hr = (Ripk3.LPS_6hr.AVG - B6.LPS_6hr.AVG),
                          LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr = (Ripk3_Casp8.LPS_6hr.AVG - Ripk3.LPS_6hr.AVG),
                          LogFC.Ripk3_Casp8.LPS_6hr.vs.B6.LPS_6hr = (Ripk3_Casp8.LPS_6hr.AVG - B6.LPS_6hr.AVG),
                          LogFC.Ripk3_Casp8.untreated.vs.Ripk3.untreated = (Ripk3_Casp8.untreated.AVG - Ripk3.untreated.AVG),
                          LogFC.Ripk3.untreated.vs.B6.untreated = (Ripk3.untreated.AVG - B6.untreated.AVG),
                          LogFC.Ripk3_Casp8.untreated.vs.B6.untreated = (Ripk3_Casp8.untreated.AVG - B6.untreated.AVG),
                          geneSymbols=geneSymbols) # why do we have a column at the end with gene Symbols??

#now look at this modified data table
#row.names(myDiffGenes.LPS) <- geneSymbols
head(myDiffGenes.LPS)
write.table(myDiffGenes.LPS, "myCluster3.1.xls", sep="\t", quote=FALSE)

#use dplyr "arrange" and "select" functions to sort by LogFC column of interest (arrange)
#and then display only the columns of interest (select) to see the most differentially expressed genes
myDiffGenes.LPS.sort.R3C8vR3 <- myDiffGenes.LPS %>%
  arrange(desc(LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr)) %>%
  select(geneSymbols, LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr, Ripk3.LPS_6hr.AVG, Ripk3_Casp8.LPS_6hr.AVG) 
head(myDiffGenes.LPS.sort.R3C8vR3)
write.table(myDiffGenes.LPS.sort.R3C8vR3, "myDiffGenes.LPS.sort.R3C8vR3.xls", sep="\t", quote=FALSE)

myDiffGenes.LPS.sort.R3vB6 <- myDiffGenes.LPS %>%
  arrange(desc(LogFC.Ripk3.LPS_6hr.vs.B6.LPS_6hr)) %>%
  select(geneSymbols, LogFC.Ripk3.LPS_6hr.vs.B6.LPS_6hr, B6.LPS_6hr.AVG, Ripk3.LPS_6hr.AVG) 
head(myDiffGenes.LPS.sort.R3vB6)
write.table(myDiffGenes.LPS.sort.R3vB6, "myDiffGenes.LPS.sort.R3vB6.xls", sep="\t", quote=FALSE)

myDiffGenes.LPS.sort.R3C8vB6 <- myDiffGenes.LPS %>%
  arrange(desc(LogFC.Ripk3_Casp8.LPS_6hr.vs.B6.LPS_6hr)) %>%
  select(geneSymbols, LogFC.Ripk3_Casp8.LPS_6hr.vs.B6.LPS_6hr, B6.LPS_6hr.AVG, Ripk3_Casp8.LPS_6hr.AVG) 
head(myDiffGenes.LPS.sort.R3C8vB6)
write.table(myDiffGenes.LPS.sort.R3C8vB6, "myDiffGenes.LPS.sort.R3C8vB6.xls", sep="\t", quote=FALSE)

#use dplyr "filter" and "select" functions to pick out genes of interest (filter)
#ways to tweek the 'select' function
#use : between two column names to select all columns between
#use 'contains', 'starts_with' or 'ends_with' to modify how you select
#can refer to columns using exact name or numerical indicator
myDiffGenes.LPS.filter <- myDiffGenes.LPS %>%
  filter(geneSymbols=="Cebpb" | geneSymbols=="Bach1") %>%
  select(geneSymbols, LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr, LogFC.Ripk3_Casp8.LPS_6hr.vs.B6.LPS_6hr) 
head(myDiffGenes.LPS.filter)

#another example using grep to match patterns - how do you filter multiple patterns of rows into one data table?
myDiffGenes.LPS.filter <- myDiffGenes.LPS %>%
  #filter(grepl('Ccl', geneSymbols)) %>%
  filter(grepl('Il', geneSymbols) | grepl('Tnf', geneSymbols) | grepl('Ccl', geneSymbols) | grepl('Cxcl', geneSymbols)) %>%
  select(geneSymbols, LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr) 
head(myDiffGenes.LPS.filter)

#first reorder and clean up the data you want to graph
#row.names(myCollapsedData.filter) <- myCollapsedData.filter[,1]
#myCollapsedData.filter <- select(myCollapsedData.filter, -geneSymbols)
#myCollapsedData.filter.transpose <- as.data.frame(t(myCollapsedData.filter))
#cellType <- c("B1a", "B1b", "Mac")
#myCollapsedData.filter.transpose <- mutate(myData.filter.transpose, cellType, 
phenotype=c("naive", "naive","naive"))
#myCollapsedData.filter.transpose

#another way to set up for graphing...don't filter before transposing
#myCollapsedData.graph <- select(myCollapsedData,-geneSymbols)
#head(myCollapsedData.graph)
#myCollapsedData.graph.transpose <- as.data.frame(t(myCollapsedData.graph))
#dim(myCollapsedData.graph.transpose)

#plot a simple bar graph - y=column of interest, x= row of interest from filtered table
ggplot(myDiffGenes.LPS.filter, aes(y=LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr)) +
  geom_bar(aes(x=geneSymbols, fill=LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr), stat="identity") +
  theme(axis.text.x=element_text(angle=-45))

#create a basic scatterplot using ggplot
ggplot(myDiffGenes.LPS, aes(x=Ripk3.LPS_6hr.AVG, y=Ripk3_Casp8.LPS_6hr.AVG)) +
  geom_point(shape=1) +
  geom_point(size=4)

# ##Volcano plots
# ##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
# gene_list$threshold = as.factor(abs(gene_list$logFC) > 2 & gene_list$P.Value < 0.05/no_of_genes)
# 
# ##Construct the volcano plot
# ggplot(data=myData, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
#   geom_point(alpha=0.4, size=1.75) +
#   opts(legend.position = "none") +
#   xlim(c(-10, 10)) + ylim(c(0, 15)) +
#   xlab("log2 fold change") + ylab("-log10 p-value")

#define a tooltip that shows gene symbol and Log2 expression data when you mouse over each data point in the plot
tooltip <- function(data, ...) {
  paste0("<b>","Symbol: ", data$geneSymbols, "</b><br>",
         "Ripk3.LPS_6hr.AVG: ", data$Ripk3.LPS_6hr.AVG, "<br>",
         "Ripk3_Casp8.LPS_6hr.AVG: ", data$Ripk3_Casp8.LPS_6hr.AVG)
}

#plot the interactive graphic
myDiffGenes.LPS %>% 
  ggvis(x= ~Ripk3.LPS_6hr.AVG, y= ~Ripk3_Casp8.LPS_6hr.AVG, key := ~geneSymbols) %>% 
  layer_points(fill = ~LogFC.Ripk3_Casp8.LPS_6hr.vs.Ripk3.LPS_6hr) %>%
  add_tooltip(tooltip)

#####################################################################################

#Step 5 HEATMAPS

#this script creates detailed heatmaps from your differentially expressed genes
#first part of script contains basic heatmap creation by reading in a file of a priori genes
#second part of script uses all the differentially expressed genes you identified previously

###################################################################################################
# generate a heatmap of differentially expressed transcripts using the entire dataset
###################################################################################################
#LPS HEATMAPS

#first let's generate a heatmap of diff expr genes in response to LPS treatment (contrast matrix.LPS - see above)
diffData.LPS.subset <- diffData.LPS[,c(7,8,9,10,11,12)]
hr <- hclust(as.dist(1-cor(t(diffData.LPS.subset), method="pearson")), method="average") #cluster rows by pearson correlation
hc <- hclust(as.dist(1-cor(diffData.LPS.subset, method="spearman")), method="complete") #cluster columns by spearman correlation
# Cut the resulting tree and create color vector for clusters.  Vary the cut height to give more or fewer clusters, or you the 'k' argument to force n number of clusters
dim(diffData.LPS.subset)
#cutree will cut to give 7 groups if you do k=7
mycl <- cutree(hr, k=4)
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9) 
mycolhc <- mycolhc[as.vector(mycl)] 
#load the gplots package for plotting the heatmap
library(gplots) 
#assign your favorite heatmap color scheme. Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]).
myheatcol <- greenred(75)
library(RColorBrewer)

#plot the hclust results as a heatmap (arg1=datatable, Rowv=what do you want in your rows?)
#heatmap.2(diffData.LPS, Rowv=as.dendrogram(hr), Colv=NA, 
          #col=myheatcol, scale="row", labRow=NA, labCol=NA,
          #density.info="none", trace="none", RowSideColors=mycolhc, 
          #cexRow=1, cexCol=1, margins=c(8,30)) 
#notice that the heatmap includes ALL the columns from your dataset

heatmap.2(diffData.LPS.subset, Rowv=as.dendrogram(hr), Colv=NA, 
          col=myheatcol, scale="row", labRow=NA, labCol=NA,
          density.info="none", trace="none", RowSideColors=mycolhc, 
          cexRow=1, cexCol=1, margins=c(8,30)) 

#now find out what your clusters are 
names(mycolhc) <- names(mycl) 
barplot(rep(10, max(mycl)), 
        col=unique(mycolhc[hr$labels[hr$order]]), 
        horiz=T, names=unique(mycl[hr$order])) # Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'.

#Now repeat this process, but with your biological replicates averaged
head(diffData.LPS)
library(limma)
colnames(diffData.LPS) <- myGroups
#cannot do this for the subset do it for whole thing first
dim(diffData.LPS)
#rownames(diffData.LPS) <- diffProbesLPS
head(diffData.LPS)
#dim(diffProbesLPS)
diffData.LPS.AVG <- avearrays(diffData.LPS, ID=colnames(diffData.LPS))
head(diffData.LPS.AVG)
dim(diffData.LPS.AVG)

#then rerun the heatmap script above using diffData.AVG as input instead of diffData

diffData.LPS.subset.AVG <- diffData.LPS.AVG[,c(4,5,6)]

hr <- hclust(as.dist(1-cor(t(diffData.LPS.subset.AVG), method="pearson")), method="complete") #cluster rows by pearson correlation
hc <- hclust(as.dist(1-cor(diffData.LPS.subset.AVG, method="spearman")), method="complete") #cluster columns by spearman correlation
# Cut the resulting tree and create color vector for clusters.  Vary the cut height to give more or fewer clusters, or you the 'k' argument to force n number of clusters
dim(diffData.LPS.subset.AVG)
#cutree will cut to give 7 groups if you do k=7
mycl <- cutree(hr, k=6)
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9) 
mycolhc <- mycolhc[as.vector(mycl)] 
#load the gplots package for plotting the heatmap
library(gplots) 
#assign your favorite heatmap color scheme. Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]).
myheatcol <- greenred(75)
library(RColorBrewer)
heatmap.2(diffData.LPS.subset.AVG, Rowv=as.dendrogram(hr), Colv=NA, 
          col=myheatcol, scale="row", labRow=NA, labCol=NA,
          density.info="none", trace="none", RowSideColors=mycolhc, 
          cexRow=1, cexCol=1, margins=c(8,30), key=T) 
#clusters 4, 3, 2, 1 from top down.
#clusters 5, 4, 3, 6, 2, 1 from top down.

#select sub-clusters of co-regulated transcripts for downstream analysis

clid <- c(3) 
ysub <- diffData.LPS.subset.AVG[names(mycl[mycl%in%clid]),] 
hrsub <- hclust(as.dist(1-cor(t(ysub), method="pearson")), method="average") 
clusterIDs <- data.frame(Labels=rev(hrsub$labels[hrsub$order]))
clusterIDs <- as.vector(t(clusterIDs))
heatmap.2(ysub, Rowv=as.dendrogram(hrsub), Colv=NA, labRow=NA, col=myheatcol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc[mycl%in%clid], margins=c(8,30)) # Create heatmap for chosen sub-cluster.

#retrieve gene symbols and entrezIDs for selected cluster and print out to an excel spreadsheet for downstream applications (i.e. GO enrichment in DAVID)
write.table(clusterIDs, "LPS.AVG.Cluster3.6groups.xls", sep="\t", quote=FALSE)

####################################################################################################################################################

#UNTREATED HEATMAPS

#first let's generate a heatmap of diff expr genes for untreated samples (contrast matrix.untreated - see above)
diffData.Untreated.subset <- diffData.Untreated[,c(1,2,3,4,5,6)]
hr <- hclust(as.dist(1-cor(t(diffData.Untreated.subset), method="pearson")), method="complete") #cluster rows by pearson correlation
hc <- hclust(as.dist(1-cor(diffData.Untreated.subset, method="spearman")), method="complete") #cluster columns by spearman correlation
# Cut the resulting tree and create color vector for clusters.  Vary the cut height to give more or fewer clusters, or you the 'k' argument to force n number of clusters
dim(diffData.Untreated.subset)
#cutree will cut to give 7 groups if you do k=7
mycl <- cutree(hr, k=6)
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9) 
mycolhc <- mycolhc[as.vector(mycl)] 
#load the gplots package for plotting the heatmap
library(gplots) 
#assign your favorite heatmap color scheme. Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]).
myheatcol <- greenred(75)
library(RColorBrewer)

#plot the hclust results as a heatmap (arg1=datatable, Rowv=what do you want in your rows?)
#heatmap.2(diffData.LPS, Rowv=as.dendrogram(hr), Colv=NA, 
#col=myheatcol, scale="row", labRow=NA, labCol=NA,
#density.info="none", trace="none", RowSideColors=mycolhc, 
#cexRow=1, cexCol=1, margins=c(8,30)) 
#notice that the heatmap includes ALL the columns from your dataset

heatmap.2(diffData.Untreated.subset, Rowv=as.dendrogram(hr), Colv=NA, 
          col=myheatcol, scale="row", labRow=NA, labCol=NA,
          density.info="none", trace="none", RowSideColors=mycolhc, 
          cexRow=1, cexCol=1, margins=c(8,30)) 
#clusters 6,2,5,4,3,1 from top down.

#now find out what your clusters are 
names(mycolhc) <- names(mycl) 
barplot(rep(10, max(mycl)), 
        col=unique(mycolhc[hr$labels[hr$order]]), 
        horiz=T, names=unique(mycl[hr$order])) # Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'.

#Now repeat this process, but with your biological replicates averaged
head(diffData.Untreated)
library(limma)
colnames(diffData.Untreated) <- myGroups
#cannot do this for the subset do it for whole thing first
dim(diffData.Untreated)
#rownames(diffData.LPS) <- diffProbesLPS
head(diffData.Untreated)
#dim(diffProbesLPS)
diffData.Untreated.AVG <- avearrays(diffData.Untreated, ID=colnames(diffData.Untreated))
head(diffData.Untreated.AVG)
dim(diffData.Untreated.AVG)

#then rerun the heatmap script above using diffData.AVG as input instead of diffData

diffData.Untreated.subset.AVG <- diffData.Untreated.AVG[,c(1,2,3)]

hr <- hclust(as.dist(1-cor(t(diffData.Untreated.subset.AVG), method="pearson")), method="complete") #cluster rows by pearson correlation
hc <- hclust(as.dist(1-cor(diffData.Untreated.subset.AVG, method="spearman")), method="complete") #cluster columns by spearman correlation
# Cut the resulting tree and create color vector for clusters.  Vary the cut height to give more or fewer clusters, or you the 'k' argument to force n number of clusters
dim(diffData.Untreated.subset.AVG)
#cutree will cut to give 7 groups if you do k=7
mycl <- cutree(hr, k=4)
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9) 
mycolhc <- mycolhc[as.vector(mycl)] 
#load the gplots package for plotting the heatmap
library(gplots) 
#assign your favorite heatmap color scheme. Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]).
myheatcol <- greenred(75)
library(RColorBrewer)
heatmap.2(diffData.Untreated.subset.AVG, Rowv=as.dendrogram(hr), Colv=NA, 
          col=myheatcol, scale="row", labRow=NA, labCol=NA,
          density.info="none", trace="none", RowSideColors=mycolhc, 
          cexRow=1, cexCol=1, margins=c(8,30), key=T) 
#clusters 4,3,2,1 from top down.

###############################################################################################
#select sub-clusters of co-regulated transcripts for downstream analysis
###############################################################################################
clid <- c(4) 
ysub <- diffData.Untreated.subset.AVG[names(mycl[mycl%in%clid]),] 
hrsub <- hclust(as.dist(1-cor(t(ysub), method="pearson")), method="average") 
clusterIDs <- data.frame(Labels=rev(hrsub$labels[hrsub$order]))
clusterIDs <- as.vector(t(clusterIDs))
heatmap.2(ysub, Rowv=as.dendrogram(hrsub), Colv=NA, labRow=NA, col=myheatcol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc[mycl%in%clid], margins=c(20,25)) # Create heatmap for chosen sub-cluster.

#retrieve gene symbols and entrezIDs for selected cluster and print out to an excel spreadsheet for downstream applications (i.e. GO enrichment in DAVID)
write.table(clusterIDs, "Untreated.AVG.Cluster4.xls", sep="\t", quote=FALSE)

###############################################################################################
#read in your own data to make a heatmap
###############################################################################################

mySelectedHeatmap <- read.delim("LPS.Cluster3.1.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, row.names=1)
mySelectedHeatmap <- mySelectedHeatmap[,c(-1,-2,-3,-4,-5,-6)]
mySelectedHeatmap <- as.matrix(mySelectedHeatmap)
#carry out hclust on the collapsed data matrix to generate a distance matrix for clustering

hr <- hclust(as.dist(1-cor(t(mySelectedHeatmap), method="pearson")), method="complete") #cluster rows by pearson correlation
hc <- hclust(as.dist(1-cor(mySelectedHeatmap, method="spearman")), method="average") #cluster columns by spearman correlation
dim(mySelectedHeatmap)
mycl <- cutree(hr, k=2)
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9) 
mycolhc <- mycolhc[as.vector(mycl)] 
#load the gplots package for plotting the heatmap
library(gplots) 
#assign your favorite heatmap color scheme. Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]).
myheatcol <- greenred(75)
library(RColorBrewer)
#Rowv=NA maintains the order of the rows you have in your data table.
heatmap.2(mySelectedHeatmap, Rowv=NA, Colv=NA, col=myheatcol, 
          scale="row", density.info="none", trace="none", labCol=NA, 
          RowSideColors=mycolhc, labRow=NULL,
          cexRow=1.5, cexCol=1, margins=c(5,30), key=T) # Creates heatmap for entire data set where the obtained clusters are indicated in the color bar.

###############################################################################################
#yet another way is to first generate a table using dplyr 'filter' and 'select' to get table for heatmap
###############################################################################################
myData.filter <- myData %>%
  filter((abs(Ecdysone.vs.PBS_18hr_gut) >= 1) | (abs(Ecdysone.vs.PBS_5hr_gut) >= 1)) %>%
  select(geneID, Ecdysone.vs.PBS_5hr_carcass, Ecdysone.vs.PBS_18hr_carcass, Ecdysone.vs.PBS_5hr_gut, Ecdysone.vs.PBS_18hr_gut)
head(myData.filter)
