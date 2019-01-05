# Introduction to this script -----------
# the goal of this script is to give you access to massive amounts of gene expression data without the need to download and analyze individual .fastq files
# this script allows you to access RNAseq data from GEO and SRA using the "all RNA-seq and ChIP-seq sample and signature search" (ARCHS4) project (https://amp.pharm.mssm.edu/archs4)
# this script also provides programmatic access to the "Library of Integrated Network-Based Cellular Signatures" (LINCS)  

# Load packages ------
library(tidyverse)
library(rhdf5)
library(edgeR)
library(slinky) #for interacting with LINCS data
library(csaw) #for CHIP-Seq analysis, but we'll use for handling SummarizedExperiment object

# load ARCHS4 database -----
# you should have already downloaded the most recent versions of mouse and human RNAseq data from ARCHS4 in hdf5 format
# you should watch this video from the ARCHS4 website to better understand this format: https://www.youtube.com/watch?v=TjkWSBQuKoE
# begin by setting the path to your archs4 data files
archs4.human <- "~/Dropbox/publicData/human_matrix.h5"
archs4.mouse <- "~/Dropbox/publicData/mouse_matrix.h5"

# 133,776 samples from human
samples.human <- h5read(archs4.human, "meta/Sample_geo_accession")
# 170,010 samples from mouse
samples.mouse <- h5read(archs4.mouse, "meta/Sample_geo_accession")

# query ARCHS4 database ----
# choose your samples based on GEO or SRA ID
samp <- c("GSM1224927","GSM1066120","GSM1224923","GSM1224929","GSM1224924","GSM1066118","GSM1066119","GSM1224925","GSM1224930","GSM1872071","GSM2282084","GSM1872064","GSM1872067","GSM1704845")
sample_locations <- which(samples.human %in% samp)
# Identify columns to be extracted from ARCHS4 database
tissue <- h5read(archs4.human, "meta/Sample_source_name_ch1")
genes <- h5read(archs4.human, "meta/genes")

# extract expression data from ARCHS4 ----
expression <- h5read(archs4.human, "data/expression", index=list(1:length(genes), sample_locations))
H5close()
rownames(expression) <- genes
colnames(expression) <- samples.human[sample_locations]
colSums(expression)
archs4.dge <- DGEList(expression)
archs4.log2.cpm <- cpm(archs4.dge, log=TRUE)
colnames(archs4.log2.cpm) <- colnames(expression)
write_tsv(archs4.log2.cpm, "expression.txt")


# Load LINCS L1000 database ----
# LINCS Phase 1 data is in GEO series GSE92742 and includes data for 1,319,138 samples
# LINCS Phase 2 data is in GEO series GSE70138 and includes data for 354,123 samples
# for each sample, the expression of 978 'landmark' genes was measured using luminex beads, which are then used to computationaly infer the expression of an additional 11,350 genes
# LINCS data is stored in a GCTx file, which is a format based on HDF5, which you can learn more about here (https://www.biorxiv.org/content/early/2017/11/30/227041)
# you can download and query these GCTx files, or query the main warehouse for LINCS data, called Clue.io

#begin by setting the path to our files for the L1000 phase 2 data (available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138)
key <- "44fa56ffbe35ee04fc291c73505058e2" #this is a user specific key that you'll need to generate for yourself by creating a free account on Clue.io
L1K.phase1.gctx <- "~/Dropbox/publicData/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"
L1K.phase1.info <- "~/Dropbox/publicData/GSE92742_Broad_LINCS_inst_info.txt"
L1K.phase2.gctx <- "~/Dropbox/publicData/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx"
L1K.phase2.info <- "~/Dropbox/publicData/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt"
#now you're ready to create the slinky object that stores all the data from the files above
slob.phase1 <- Slinky(key, L1K.phase1.gctx, L1K.phase1.info)
slob.phase2 <- Slinky(key, L1K.phase2.gctx, L1K.phase2.info)

# Query L1000 Phase 1 data by accessing clue.io ----
# since phase 1 contains so much data, it can be hard to parse these files locally on your laptop
# instead, we'll query the Clue.io API
# First, explore the Phase 1 data based on metadata in the info.txt file
phase1meta.df <- as.tibble(metadata(slob.phase1))
glimpse(phase1meta.df)

#take a look at the data grouped by perturbagen
purturbagen.phase1 <- phase1meta.df %>%
  group_by(pert_iname) %>%
  summarize(n())

#take a look at the data grouped by cell line
cellLine.phase1 <- phase1meta.df %>%
  group_by(cell_id) %>%
  summarize(n())

#take a look at the data grouped by timepoint
timepoint.phase1 <- phase1meta.df %>%
  group_by(pert_time) %>%
  summarize(n())

#use the 'loadL1K' function from Slinky to create an summarizedExperiment object directly from the Clue.io query
phase1 <- loadL1K(slob.phase1, where_clause = list(pert_type = "trt_cp", 
                                                   pert_iname = "quinpirole", 
                                                   cell_id = "A375", 
                                                   is_gold = TRUE), inferred = FALSE)

#convert the summarizedExperiment into a DGEList object using the 'asDGEList' function from the csaw package
phase1.DGE <- asDGEList(phase1)

# Query L1000 Phase 2 data directly from local files ----
# this approach should be used for the Phase 2 data, which is not yet on Clue.io
# begin by exploring Phase 2 based on metadata in the info.txt file
phase2meta.df <- as.tibble(metadata(slob.phase2))
glimpse(phase2meta.df)

#take a look at the data grouped by perturbagen
purturbagen.phase2 <- phase2meta.df %>%
  group_by(pert_iname) %>%
  summarize(n())

#take a look at the data grouped by cell line
cellLine.phase2 <- phase2meta.df %>%
  group_by(cell_id) %>%
  summarize(n())

#take a look at the data grouped by timepoint
timepoint.phase2 <- phase2meta.df %>%
  group_by(pert_time) %>%
  summarize(n())

# Identify samples of interest based on metadata file
col.ids <- which(metadata(slob.phase2)$pert_type == "trt_cp" & 
                  metadata(slob.phase2)$pert_iname == "aminoguanidine" & 
                  metadata(slob.phase2)$cell_id == "A375")

#retrieve data for these samples from the .gctx file
L1k.data <- readGCTX(slob.phase2[, col.ids])
L1k.dge <- DGEList(L1k.data)
L1k.log2.cpm <- cpm(L1k.dge)


