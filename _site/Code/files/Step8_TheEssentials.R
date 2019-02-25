# Introduction to this script ----
# This script brings together all the condensed 'essential' parts of the scripts you've been working with during the course

# Reproducibility ----
library(containerit)
dockerfile_object <- dockerfile()

# Step 1 ----
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(biomaRt) # provides access to a wealth of annotation info
targets <- read_tsv("../../studyDesign.txt")# read in your study design
path <- file.path("../readMapping", targets$sample, "abundance.h5") # set file paths to your mapped data
targets <- mutate(targets, path) # add paths to your study design (only necessary for Sleuth)
Hs.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl") # select 'mart' from biomaRt for annotations
Tx <- getBM(attributes=c('ensembl_transcript_id_version', # get gene symbols for each transcript ID
                         'external_gene_name'),
            mart = Hs.anno)
Tx <- as_tibble(Tx) # convert this annotation mapping file to a tibble (the tidyverse version of a dataframe)
Tx <- dplyr::rename(Tx, target_id = ensembl_transcript_id_version, gene_name = external_gene_name) # change some names
Txi_gene <- tximport(path, #reading kallisto data into R
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, 
                     countsFromAbundance = "lengthScaledTPM")
myCPM <- as_tibble(Txi_gene$abundance, rownames = "geneSymbol") # these are you counts after adjusting for transcript length
myCounts <- as_tibble(Txi_gene$counts, rownames = "geneSymbol") # these are your transcript per million (TPM) values, or counts per million (CPM) if you collapsed data to gene level

dockerfile_object
write(dockerfile_object, file = ".dockerfile")

