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
# for this part of the class you'll use your abundance data ('myTPM') that we generated in the last class
# make sure you have this object already in your work environment
# if you don't, go back to the Step2 script and generate it
head(myTPM)
colnames(myTPM)
#give your data some more informative column names
colnames(myTPM) <- labels
#convert the entire datamatrix to Log2 scale
myTPM <- log2(myTPM + 1)
# Now we need to convert our datamatrix to a dataframe
class(myTPM)
myTPM.df <- as.data.frame(myTPM)

# use dplyr 'mutate' function to add new columns based on existing data -------
myTPM.df <- mutate(myTPM.df,
                   untreated.AVG = (control.1 + control.2 + control.3)/3, 
                   wt_crypto.AVG = (wt_crypto.1 + wt_crypto.2 + wt_crypto.3)/3,
                   trans_crypto.AVG = (trans_crypto.1 + trans_crypto.2 + trans_crypto.3)/3,
                   #now make columns comparing each of the averages above that you're interested in
                   LogFC.wt_crypto_vs_untreated = (wt_crypto.AVG - untreated.AVG),
                   LogFC.trans_crypto_vs_untreated = (trans_crypto.AVG - untreated.AVG),
                   geneSymbol = row.names(myTPM.df))
#why is this type of approach to managing a spreadsheet useful?

#now look at this modified data table
head(myTPM.df)
myTPM.df.subset <- myTPM.df[,c(15,10:14)]

# using dplyr 'arrange' and 'select' to sort your dataframe based on any variable -----
# first, we'll use dplyr "arrange" function to sort rows based on the values in a column of interest
# then we'll display 'select' only the columns we're interested in seeing
myTPM.sort <- myTPM.df %>%
  dplyr::arrange(desc(LogFC.wt_crypto_vs_untreated)) %>% #note that this is the first time you've seen the 'pipe' operator
  dplyr::select(geneSymbol, LogFC.wt_crypto_vs_untreated)

head(myTPM.sort)

# use dplyr "filter" and "select" functions to pick out genes of interest  ----
#ways to tweek the 'select' function
#use : between two column names to select all columns between
#use 'contains', 'starts_with' or 'ends_with' to modify how you select
#can refer to columns using exact name or numerical indicator
#use boolean operators such as '&' (and), '|' (or), '==' (equal to), '!' (not)
myTPM.filter <- myTPM.df %>%
  dplyr::filter(geneSymbol=="IL1A" | geneSymbol=="CXCL10" | geneSymbol=="RIPK3" | geneSymbol=="CASP8") %>%
  dplyr::select(geneSymbol, untreated.AVG, wt_crypto.AVG)
myTPM.filter

# you can also filter based on any regular expression
myTPM.grep <- myTPM.df %>%
  dplyr::filter(grepl('CXCL|CXCR', geneSymbol)) %>%
  dplyr::select(geneSymbol, untreated.AVG, wt_crypto.AVG)

myTPM.grep

# make an interactive table ----
library(DT)
datatable(myTPM.df.subset, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'my cool table)',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(1:6), digits=3)

# make a simple scatter plot -----
ggplot(myTPM.df, aes(x=untreated.AVG, y=wt_crypto.AVG)) +
  geom_point(shape=3) +
  geom_point(size=1)

# make an interactive scatter plot ----
tooltip <- function(data, ...) {
  paste0("<b>","Symbol: ", data$geneSymbol, "</b><br>",
         "untreated.AVG: ", data$untreated.AVG, "<br>",
         "wt_crypto.AVG: ", data$wt_crypto.AVG)
}

#plot the interactive graphic
myTPM.df %>% 
  ggvis(x= ~untreated.AVG, y= ~wt_crypto.AVG, key := ~geneSymbol) %>% 
  #layer_points(fill = ~LogFC.B6.LPS_vs_B6.untreated) %>%
  add_tooltip(tooltip)
#notice that your R session is always active during an interactive session
