# Introduction to this script -----------
#This script carries out the steps involved in analysis of RNAseq data.  
#depending on your data and interests, different parts of this script may not apply to you

# Loading packages -----------------------------------------------
#begin by loading the packages required for RNAseq data
library(Rsubread) #package that I use for read alignment and summarization. NOTE: this package ONLY available for Mac
library(limma) #comprehensive package for using linear models to analyze gene expression data
library(edgeR) #comprehensive package for analysis of RNAseq data
library(ShortRead) #working with short sequences
library(matrixStats) #stats and calculations
library(ggplot2) #great package for publication quality graphs
library(org.Mm.eg.db) #our organism-specific database package that we'll use for annotation
library(AnnotationDbi) #

# Study design --------------------------------------------
#In this section, you will layout the design of your experiment using a text file
#begin by reading in a text file where columns contain metadata and rows represent samples
targets <- read.delim("Beiting_studyDesign.txt", row.names=NULL)
#capture information about treatment groups from this file
groups <- paste(targets$genotype, targets$treatment, sep=".")
#turn this variable into a 'factor'
groups <- factor(groups)
#create some more human-readable labels for your samples using the info in your studyDesign file
sampleLabels <- paste(targets$genotype, targets$treatment, targets$rep, sep=".")

#set-up your experimental design
design <- model.matrix(~0+groups)
colnames(design) <- levels(groups)
design

# QC metrics for fastq ------------------------------------------
# #you can check read quality using shortRead package
# #but I usually find it is better to do this on the sequencer or Illumina's BaseSpace website
# myFastq <- targets$fastq
# #collecting statistics over the files
# qaSummary <- qa("B6-WT-untreat-rep2_S2_mergedLanes_read1.fastq", type="fastq")
# #create and view a report
# browseURL(report(qaSummary))

# Build Index from genome fasta --------------------------------
# expect this to take about 20 min on 8G RAM for mouse or genome
#you must have already downloaded the fasta file for your genome of interest and have it in the working directory
#this only needs to be done once, then index can be reused for future alignments
buildindex(basename="mouse",reference="Mus_musculus.GRCm38.dna.primary_assembly.fa")

# Align reads -----------------------------------------
#align your reads (in the fastq files) to your indexed reference genome that you created above
#expect this to take about 45min for a single fastq file containing 25 million reads
#the output from this is a .bam file for each of your original fastq files
reads1 <- targets$fastq[9]
reads2 <- targets$fastq[21] 
align(index="mouse", readfile1=reads1, readfile2=reads2, input_format="gzFASTQ",output_format="BAM",
      output_file="alignmentResultsPE_sample9.BAM", tieBreakHamming=TRUE,unique=TRUE,indels=5, nthreads=8)

# Summarize reads -------------------------------
#use the 'featureCounts' function to summarize read counts to genomic features (exons, genes, etc)
#will take about 1-2min per .bam file.
#for total transcriptome data summarized to mouse/human .gtf, expect about 50-60% of reads summarize to genes (rest is non-coding)
#read in text file with bam file names
MyBAM <- read.delim("BAMfiles_names.txt", header=T)
MyBAM <- as.character(MyBAM[,1])
#summarize aligned reads to genomic features (i.e. exons)
fc <- featureCounts(files=MyBAM, annot.ext="Mus_musculus.GRCm38.79.gtf", isGTFAnnotationFile=TRUE, GTF.featureType = "exon",
                    GTF.attrType="gene_id", useMetaFeatures=TRUE, isPairedEnd=TRUE, requireBothEndsMapped=TRUE, strandSpecific=2, nthreads=8)
#use the 'DGEList' function from EdgeR to make a 'digital gene expression list' object
DGEList <- DGEList(counts=fc$counts, genes=fc$annotation)
DGEList

# Looking at your data as a DGEList object -------------------------
#If you provided raw .fastq files at the start of the class, then you will start you analysis here
#by reading in a 'DGEList' object (created above by the edgeR package)
#we created this DGEList object for you
#begin by loading the DGEList object
load("DGEList") #once loaded, take a look at this contents of this object

#retrieve all your gene/transcript identifiers from this DGEList object
myEnsemblIDs <- DGEList$genes$GeneID

#this next bit is not essential, but it's useful to graph your data before and after Log2 transformation
counts <- DGEList$counts
counts.Log2 <- log2(counts + 0.5)
countsSD <- transform(counts, SD=rowSds(counts, na.rm=TRUE), AVG=rowMeans(counts), MED=rowMedians(counts))
countsSD.Log2 <- transform(counts.Log2, SD=rowSds(counts.Log2, na.rm=TRUE), AVG=rowMeans(counts.Log2), MED=rowMedians(counts.Log2))

#now graph this data.
ggplot(countsSD, aes(x=SD, y=MED)) +
  geom_point(shape=1) +
  geom_point(size=4)
#the graph shows the non-Log2 data (median vs StdDev for each row) 
#notice how the highly expressed genes have highest standard deviation.
#why might this be problem? 

#now graph the Log2 data, to see how this basic transformation changes the structure of your data 
ggplot(countsSD.Log2, aes(x=SD, y=MED)) +
  geom_point(shape=1) +
  geom_point(size=4)
#notice that the lowest expressed genes now have the highest StdDev

# Normalize using mean-variance relationship ---------------
normData.unfiltered <- voom(DGEList, design, plot=TRUE)
exprs.unfiltered <- normData.unfiltered$E
exprs.matrix.unfiltered <- as.matrix(exprs.unfiltered)
#note that because you're now working with Log2 CPM, much of your data will be negative number (log2 of number smaller than 1 is negative) 
head(exprs.matrix.unfiltered)

#if you need RPKM for your unfiltered, they can generated as follows
#Although RPKM are commonly used, not really necessary since you don't care to compare two different genes within a sample
rpkm.unfiltered <- rpkm(DGEList, DGEList$genes$Length)
rpkm.unfiltered <- log2(rpkm.unfiltered + 0.5)

# Filtering your dataset ----------------
#Only keep in the analysis those genes which had >10 reads per million mapped reads in at least two samples
filter.test <- rowSums(cpm(DGEList) > 10) >= 2
DGEList.filtered <- DGEList[filter.test,]
dim(DGEList.filtered)
rpkm.filtered <- rpkm(DGEList.filtered, DGEList.filtered$genes$Length) #if you prefer, can use 'cpm' instead of 'rpkm' here
rpkm.filtered <- log2(rpkm.filtered + 1)

dim(exprs.matrix.unfiltered)

#Use Voom again to normalize this filtered dataset
normData.filtered <- voom(DGEList.filtered, design, plot=TRUE)
exprs <- normData.filtered$E
exprs.matrix.filtered <- as.matrix(exprs)
head(exprs.matrix.filtered)

# Annotating your data -----------------
#If we want to know what kinds of data are retriveable via the 'select' command, look at the columns of the annotation database
columns(org.Mm.eg.db)
#If we want to know what kinds of fields we could potentially use as keys to query the database, use the 'keytypes' command
keytypes(org.Mm.eg.db)
#transform you identifiers to entrezIDs
myAnnot.unfiltered <- AnnotationDbi::select(org.Mm.eg.db, keys=rownames(exprs.matrix.unfiltered), keytype="ENSEMBL", columns=c("ENTREZID", "GENENAME", "SYMBOL"))
myAnnot.filtered <- AnnotationDbi::select(org.Mm.eg.db, keys=rownames(exprs.matrix.filtered), keytype="ENSEMBL", columns=c("ENTREZID", "GENENAME", "SYMBOL"))
resultTable.unfiltered <- merge(myAnnot.unfiltered, exprs.matrix.unfiltered, by.x="ENSEMBL", by.y=0)
resultTable.filtered <- merge(myAnnot.filtered, exprs.matrix.filtered, by.x="ENSEMBL", by.y=0)
head(resultTable.unfiltered)
#add more appropriate sample names as column headers
colnames(resultTable.unfiltered) <- c("Ensembl", "entrez", "name", "symbol", sampleLabels)
colnames(resultTable.filtered) <- c("Ensembl", "entrez", "name", "symbol", sampleLabels)
#now write these annotated datasets out
write.table(resultTable.unfiltered, "normalizedUnfiltered.txt", sep="\t", quote=FALSE)
write.table(resultTable.filtered, "normalizedFiltered.txt", sep="\t", quote=FALSE)
head(resultTable.filtered)
#congrats! You now have normalized and annotated RNAseq data
#you also have several ways to modify this data (Log2, CPM, RPKM, filtering)

# OPTIONAL: Starting with RPKM or CPM tables -----------------
#There may be times where you only have access to processed data, not the original fastq or bam files
#in this case you can use text files with RPKM or CPM data to begin your analysis
#reading in the RPKM or CPM data from a simple text file
myData <- read.delim("myData.txt", header=T)
#take a look at this data
head(myData)

#create new objects in R that will capture the gene symbols and IDs
mySymbols <- myData$Symbol
myGeneID <- myData$Gene

#make another version of the data that is only numeric gene expression data, by removing all columns that have text (gene symbols, IDs, descriptions, etc)
myData.matrix <- myData[,-1:-2]
myData.matrix <- as.matrix(myData.matrix)
head(myData.matrix) #notice in your environment how this object is different than the 'myData' object you initially created

#now you can filter to remove lowly expressed transcripts
filter.test <- rowSums((myData.matrix) > 10) >= 2
myData.matrix.filtered <- myData.matrix[filter.test,]
dim(myData.matrix.filtered)
#you don't even need to think about annotation, since if someone gave you a table of RPKM data, they had probably already annotated the data

# OPTIONAL: alternative annotation approach for non-model organisms -----------
# #If there isn't an R database package for your organism,
# #you can still annotate your data using the biomaRt package
# library(biomaRt)
# #list all the available 'marts' to choose from
# listMarts()
# #select your mart of interest
# vectorBase <- useMart("vb_gene_mart_1502")
# #list the databases available in that mart
# listDatasets(vectorBase)
# #select the database you want to work with
# aedesDB <- useDataset("aaegypti_eg_gene", mart=vectorBase)
# #list the various filters, or keys, that can used to access this database
# listFilters(aedesDB)
# #list all the attributes that can be retrieved from the database
# listAttributes(aedesDB)
# #pull out your rownames of your dataset to use as one of the filters
# myEnsemblIDs.filtered <- rownames(rpkm.filtered)
# myEnsemblIDs.unfiltered <- rownames(rpkm.unfiltered)
# #attach annotation information to your file
# myAnnot.filtered <- getBM(attributes = c("external_gene_id", "ensembl_gene_id", "description"),
#                           filters = "ensembl_gene_id",   #the kind of data you'll use as keys
#                           values = myEnsemblIDs.filtered,   #the exact data you'll use as keys
#                           mart=aedesDB)
# myAnnot.unfiltered <- getBM(attributes = c("external_gene_id", "ensembl_gene_id", "description"),
#                             filters = "ensembl_gene_id",
#                             values = myEnsemblIDs.unfiltered,
#                             mart=aedesDB)
# #merge the annotation information with the expression data
# resultTable.filtered <- merge(myAnnot.filtered, rpkm.filtered, by.x="ensembl_gene_id", by.y=0)
# resultTable.unfiltered <- merge(myAnnot.unfiltered, rpkm.unfiltered, by.x="ensembl_gene_id", by.y=0)

