# Introduction to this script -----------
# this script walks thorough some basic data wrangling for dealing with any dataframe
# we'll also continue to create publication quality graphics
# we'll start the script with your abundance data from the 

# Load packages ------
library(ggplot2) #graphing package that employs a 'grammar of graphics' approach
library(reshape2) #data manipulation package
library(dplyr) #data manipulation package
library(ggvis) #for making interactive (dynamic) graphs in R, also follows the grammar of graphics philosophy
options(digits=3)

# Data -------
# for this part of the class you'll use your normalized and filtered data in log2 cpm
# make sure you have this object already in your work environment
# if you don't, go back to the Step2 script and generate it
head(log2.cpm.filtered.norm)
#give your data some more informative column names
colnames(log2.cpm.filtered.norm) <- sampleLabels
# Now we need to convert our datamatrix to a dataframe
class(log2.cpm.filtered.norm)
mydata.df <- as.data.frame(log2.cpm.filtered.norm)

# use dplyr 'mutate' function to add new columns based on existing data -------
mydata.df <- mutate(mydata.df,
                   untreated.AVG = (control_rep1 + control_rep2 + control_rep1)/3, 
                   wt_crypto.AVG = (wt_rep1 + wt_rep2 + wt_rep3)/3,
                   trans_crypto.AVG = (transgenic_rep1 + transgenic_rep2 + transgenic_rep3)/3,
                   #now make columns comparing each of the averages above that you're interested in
                   LogFC.wt_crypto_vs_untreated = (wt_crypto.AVG - untreated.AVG),
                   LogFC.trans_crypto_vs_untreated = (trans_crypto.AVG - untreated.AVG),
                   geneSymbol = row.names(mydata.df))
#why is this type of approach to managing a spreadsheet useful?

#now look at this modified data table
head(mydata.df)

# using dplyr 'arrange' and 'select' to sort your dataframe based on any variable -----
# first, we'll use dplyr "arrange" function to sort rows based on the values in a column of interest
# then we'll display 'select' only the columns we're interested in seeing
mydata.sort <- mydata.df %>%
  dplyr::arrange(desc(LogFC.wt_crypto_vs_untreated)) %>% #note that this is the first time you've seen the 'pipe' operator
  dplyr::select(geneSymbol, LogFC.wt_crypto_vs_untreated)


# use dplyr "filter" and "select" functions to pick out genes of interest  ----
#ways to tweek the 'select' function
#use : between two column names to select all columns between
#use 'contains', 'starts_with' or 'ends_with' to modify how you select
#can refer to columns using exact name or numerical indicator
#use boolean operators such as '&' (and), '|' (or), '==' (equal to), '!' (not)
mydata.filter <- mydata.df %>%
  dplyr::filter(geneSymbol=="MX1" | geneSymbol=="IRF1" | geneSymbol=="OAS2" | geneSymbol=="IRF3") %>%
  dplyr::select(geneSymbol, untreated.AVG, wt_crypto.AVG)
mydata.filter

# you can also filter based on any regular expression
mydata.grep <- mydata.df %>%
  dplyr::filter(grepl('CXCL|IFI', geneSymbol)) %>%
  dplyr::select(geneSymbol, untreated.AVG, wt_crypto.AVG)

mydata.grep

# make an interactive table ----
library(DT)
datatable(mydata.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'my cool table)',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(1:6), digits=3)

# make a simple scatter plot -----
ggplot(mydata.df, aes(x=untreated.AVG, y=wt_crypto.AVG)) +
  geom_point(shape=3) +
  geom_point(size=1)

# make an interactive scatter plot ----
tooltip <- function(data, ...) {
  paste0("<b>","Symbol: ", data$geneSymbol, "</b><br>",
         "untreated.AVG: ", data$untreated.AVG, "<br>",
         "wt_crypto.AVG: ", data$wt_crypto.AVG)
}

#plot the interactive graphic
mydata.df %>% 
  ggvis(x= ~untreated.AVG, y= ~wt_crypto.AVG, 
        key := ~geneSymbol,
        size := input_slider(10, 100),
        opacity := input_slider(0, 1)) %>% 
        #layer_points(fill := ~LogFC.wt_crypto_vs_untreated)) %>%
  add_tooltip(tooltip)
#notice that your R session is always active during an interactive session
