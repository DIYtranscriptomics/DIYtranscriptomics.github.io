# load libraries ----
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(biomaRt) # an alternative for annotation

# QUESTION 1 ----
# retrieve the requested annotation data for the ferret
listMarts() #default host is ensembl.org, and most current release of mammalian genomes
#choose the 'mart' you want to work with
myMart <- useMart(biomart="ENSEMBL_MART_ENSEMBL")
#take a look at all available datasets within the selected mart
available.datasets <- listDatasets(myMart)
#now grab the ensembl annotations for ferret
ferret.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "mpfuro_gene_ensembl")
ferret.attributes <- listAttributes(ferret.anno)

Tx.ferret <- getBM(attributes=c('ensembl_transcript_id',
                                'start_position',
                                'end_position',
                                'external_gene_name',
                                'description',
                                'entrezgene_id',
                                'pfam'),
                   mart = ferret.anno)

Tx.ferret <- as_tibble(Tx.ferret)

# QUESTION #2 ----
# retrieve promoter sequences for selected genes

# check out the help documentation for the getSequence function.  
?getSequence #note that the example code in the help doc is executable!

#get your ferret data from BiomaRt
ferret.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "mpfuro_gene_ensembl")

#get your promoter sequences from this annotation data
seq <- getSequence(id = c("MX1", "IFIT2", "OAS2", "IRF1", "IFNAR1", "MX1"),
                  type = "hgnc_symbol",
                  seqType = "gene_flank",
                  upstream = 1000,
                  mart = ferret.anno)

