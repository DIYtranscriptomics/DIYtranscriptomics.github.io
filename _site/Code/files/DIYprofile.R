#' Extract counts per million (CPM) from a DGEList object, convert to a tibble, pivot_longer to make tidy, and then plot with ggplot to see distribution of CPM
#'
#' @param data A object of class DGEList.
#' @param samples A character vector of sample names.
#' @return A violin plot of CPM (y axis) for each sample (x axis)
#' @export
DIYprofile <- function(data, samples) {
  log2.cpm.filtered <- edgeR::cpm(data, log=TRUE)
  log2.cpm.filtered.df <- tibble::as_tibble(log2.cpm.filtered, rownames = "geneID")
  colnames(log2.cpm.filtered.df) <- c("geneID", samples)
  log2.cpm.filtered.df.pivot <- tidyr::pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                                    cols = -1, # column names to be stored as a SINGLE variable
                                                    names_to = "samples", # name of that new variable (column)
                                                    values_to = "expression") # name of new variable (column) storing all the values (data)

  ggplot2::ggplot(log2.cpm.filtered.df.pivot) +
    aes(x=samples, y=expression, fill=samples) +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    stat_summary(fun = "median",
                 geom = "point",
                 shape = 95,
                 size = 10,
                 color = "black",
                 show.legend = FALSE) +
    labs(y="log2 expression", x = "sample",
         title="Log2 Counts per Million (CPM)",
         subtitle="filtered, non-normalized",
         caption=paste0("produced on ", Sys.time())) +
    theme_bw()

}
