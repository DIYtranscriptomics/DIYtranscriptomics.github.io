#' Gene expression data for human cutaneous leishmaniasis
#'
#' A DGEList object containing counts of reads mapping (via Kallisto) to human genes
#'
#' @format A data frame with 35643 rows (genes) and 10 variables (samples):
#' \describe{
#'   \item{HS..}{healthy subject}
#'   \item{CL..}{cutaneous leishmaniasis}
#'   ...
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/bioproject/PRJNA525604}
"myDGEList"

#' Sample metadata for each patient sample
#'
#' @format A tibble containing 10 observations of 3 variables:
#' \describe{
#'   \item{sample}{sample ID}
#'   \item{sra_accession}{accession number for the sequence read archive}
#'   \item{group}{healthy or disease (cutaneous leishmaniasis)}
#'   ...
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/bioproject/PRJNA525604}
"targets"

