#Begin by loading the Oligo package -- our main tool for analyzing Affy array data
library(oligo)

#read in the raw .CEL files (will take ~5-10min)
rawData <- read.celfiles(list.celfiles())
rawData

#normalize using RMA (will take <5min)
normData <- rma(rawData, target="core") #you have the option of seeting target="full" or "extended" if working with Exon arrays
normData

#choose color scheme for graphs
cols <- topo.colors (n=84, alpha=1)

## make different plots of non-normalized data: histogram and boxplot are most useful: others include density, pairs, MAplot, sampleRelation
sampleLabels <- read.delim("Bale_sampleLabels.txt", sep="\t", fill=TRUE, header=FALSE, as.is=TRUE)
sampleLabels <- as.character(t(sampleLabels))
hist(rawData, xlab = "log2 expression", main = "non-normalized data - histograms", col=cols)
boxplot(rawData, names = sampleLabels, ylab = "non-normalized log2 expression", main = "non-normalized data - boxplots", col=cols)
#pairs(rawData)
#MAplot(rawData)

###############################################################################################
#filter out control probes
################################################################################################
#all of Affy's Gene ST arrays have many control probes that can make output files confusing
#One option is to remove these first using the 'getMainProbes' function from the affycoretools package
library(affycoretools)
normData.main <- getMainProbes(normData)
#now extract the expression data from the eset filtered eset object you just created
normData.main.matrix <- exprs(normData.main)
dim(normData.main.matrix

#an additional filtering step can be carried out using the 'nsFilter' function 
#allows you to filter out low variance data, as well as the control probes mentioned above
library(genefilter)
#first, swap the annotation slot in the normData object
annotation(normData) <- "mogene20sttranscriptcluster.db"
filtered_geneList <- nsFilter(normData, var.func=IQR, var.filter=TRUE, var.cutoff=0.5, filterByQuantile=TRUE, feature.exclude="^AFFX")

head(filtered_geneList)
# extract the ExpressionSet from this filtered list.  If you started with an expression set object, then you would pick up at this point and skip the nsFilter step
filtered.eset <- filtered_geneList$eset
head(filtered.eset)
dim(filtered.eset)

#now convert to a datamatrix that will contain only the probes after filtering
filtered.matrix <- as.matrix(filtered.eset)
probeList <- rownames(filtered.matrix)

###############################################################################################
#Annotate your data
################################################################################################
library(annotate)
#the pd.mogene.2.0.st maps probes to probesets during the normalization/summarization step
library(pd.mogene.2.0.st)
#the mogene20sttranscriptcluster.db package maps probesets to genes in the annotation step
library(mogene20sttranscriptcluster.db)
#If we want to know what kinds of data are retriveable via the 'select' command, look at the columns of the annotation database
columns(mogene20sttranscriptcluster.db)
#If we want to know what kinds of fields we could potentially use as keys to query the database, use the 'keytypes' command
keytypes(mogene20sttranscriptcluster.db)

#first, annotate you unfiltered data
myAnnot.main <- select(mogene20sttranscriptcluster.db, keys=rownames(normData.main.matrix), keytype="PROBEID", columns=c("ENTREZID", "SYMBOL", "GENENAME"))
resultTable.main <- merge(myAnnot.main, normData.main.matrix, by.x="PROBEID", by.y=0)
head(resultTable.main)
#now annotate the variance filtered data
myAnnot <- select(mogene20sttranscriptcluster.db, keys=rownames(filtered.matrix), keytype="PROBEID", columns=c("ENTREZID", "SYMBOL", "GENENAME"))
resultTable.filtered <- merge(myAnnot, filtered.matrix, by.x="PROBEID", by.y=0)
head(resultTable.filtered)
#add more appropriate sample names as column headers
colnames(resultTable.main) <- sampleLabels
colnames(resultTable.filtered) <- sampleLabels

#now write these annotated datasets out
write.table(resultTable.main, "normalizedUnfiltered.txt", sep="\t", quote=FALSE)
write.table(resultTable.filtered, "normalizedFiltered.txt", sep="\t", quote=FALSE)
