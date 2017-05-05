# for this script to work, you will need to output two files from Illumima's GenomeStudio software 
# the first file should contain non-normalized, non-background subtracted probe-level data with columns for at least the following: Avg_Signal, Avg_NBeads, BEAD_STD and Detection Pval
# the second file should contain probe set data for controls

# being by loading the lumi package
# if using Mouse arrays, change 'Mouse' to 'Mouse' for the mapping and .db packages below
library(lumi)
library(lumiMouseIDMapping)
library(lumiMouseAll.db)

###############################################################################################
# Read in the GenomeStudio output files using the 'LumiR' and 'addControlData2lumi' functions
# This will create 'LumiBatch' objects from the control and experimental probe data
###############################################################################################
rawData <- lumiR("FinalReport_samples_probes_NoNorm_NoBkrnd.txt", convertNuID = TRUE, sep = NULL, detectionTh = 0.01, na.rm = TRUE, lib = "lumiMouseIDMapping")
rawData

# Read control probe data into a separate LumiBatch and take a look at these controls
rawData <- addControlData2lumi("FinalReport_controls_probes_NoNorm_NoBkrnd.txt", rawData)
rawData
controlData <- getControlData(rawData)
getControlType(controlData)
plotStringencyGene(rawData, lib = NULL, slideIndex = NULL, addLegend = TRUE, logMode = TRUE)
plotControlData(rawData, type = 'HOUSEKEEPING', slideIndex = NULL, logMode = TRUE, new = TRUE)
plotHousekeepingGene(rawData)

###############################################################################################
# QC check of data
###############################################################################################
summary(rawData, 'QC')

# choose color scheme for graphs
cols.ALL <- topo.colors(n=16, alpha=1)

hist(rawData, xlab = "log2 expression", main = "non-normalized data - histograms")
boxplot(rawData, ylab = "non-normalized log2 expression", main = "non-normalized data - boxplots", col=cols.ALL)
# pairs(rawData)
# MAplot(rawData)


###############################################################################################
# carry out all preprocessing steps (lumiQ, lumiB, lumiN, lumiT) using the 'lumiExpresso' function which encapsulates all these preprocessing steps
###############################################################################################
# options for normalization in Lumi include: "loess', 'quantile', 'vsn' and the default, which is robust spline normalization ('rsn') which combines the features of quantile and loess
normData <- lumiExpresso(rawData, QC.evaluation=TRUE, normalize.param=list(method='rsn'))


###############################################################################################
# repeat the summarization of QC data and graphs on the normalized data
###############################################################################################
summary(normData, 'QC')
hist(normData, xlab = "log2 expression", main = "non-normalized data - histograms", col=cols.ALL)
boxplot(normData, ylab = "non-normalized log2 expression", main = "non-normalized data - boxplots", col=cols.ALL)

######################################################
# filter out probes that don't change (low variance)
# remove duplicate genes based on entrez ID
# remove coding sequences with no entrezIDs
######################################################
library(genefilter)
filtered_geneList <- nsFilter(normData, require.entrez=TRUE, remove.dupEntrez=TRUE, var.func=IQR, var.filter=TRUE, var.cutoff=0.5, filterByQuantile=TRUE)
head(filtered_geneList)
# extract the ExpressionSet from this filtered list
filtered.eset <- filtered_geneList$eset
head(filtered.eset)

# now convert to a datamatrix that will contain only the probes after filtering
filtered.matrix <- as.matrix(filtered.eset)
probeList <- rownames(filtered.matrix)
head(filtered.matrix)


###############################################################################################
# output this filtered data to a table that can be used in Excel or other programs
# This file contains genes after normalization and filtering
###############################################################################################
# first need to get all the gene symbols and entrez IDs so you can put everything together
library(annotate)
keytypes(lumiMouseAll.db)
myGenesAll <- getSYMBOL(probeList, "lumiMouseAll.db")
myGenesAll <- as.matrix(myGenesAll)
myEntrezAll <- getEG(probeList, "lumiMouseAll.db")
myEntrezAll <- as.matrix(myEntrezAll)
write.table(cbind(myGenesAll, myEntrezAll, filtered.matrix),"normalizedFilteredData.txt", sep="\t", quote=FALSE)

###############################################################################################
# also output unfiltered expression data.  
###############################################################################################
write.exprs(normData, file='robustSplineNorm.txt')
normData.delim <- read.delim("robustSplineNorm.txt", header=TRUE, row.names=1)
normData.matrix <- as.matrix(normData.delim)
probeList2 <- rownames(normData.matrix)
library(annotate)
symbols <- getSYMBOL(probeList2, "lumiMouseAll.db")
entrezIDs <- getEG(probeList2, "lumiMouseAll.db")
write.table(cbind(symbols, entrezIDs, normData.matrix), "unfiltered_expressionData.txt", sep="\t", quote=FALSE)

