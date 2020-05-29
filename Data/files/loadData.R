# Introduction ----
# The goal of this short script to it make sure you are able to read the count data and the study design for the COVID hackdash into your R environment

# We only need a single package
library(tidyverse)
# Begin by reading in study design that includes all info for both human and ferret samples
targets <- read_tsv("covid_metadata.txt")
# then read in the human covid data and convert to a matrix with gene symbols as rownames
human_covid_data <- read_tsv("GSE147507_RawReadCounts_Human.tsv")
human_covid_data <- as.matrix(column_to_rownames(human_covid_data, "X1"))
# repeat for the ferret covid data
ferret_covid_data <- read_tsv("GSE147507_RawReadCounts_Ferret.tsv")
ferret_covid_data <- as.matrix(column_to_rownames(ferret_covid_data, "X1"))

# Now proceed with your exploration and analysis of the data!