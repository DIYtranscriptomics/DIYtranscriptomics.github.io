# Introduction to this script -----------
# this script walks thorough some basic data wrangling for dealing with any dataframe
# we'll also continue to create publication quality graphics
# we'll start the script with your abundance data from the 

# Load packages ------
library(tidyverse)
library(ggplot2) #graphing package that employs a 'grammar of graphics' approach
library(reshape2) #data manipulation package
library(dplyr) #data manipulation package
library(ggvis) #for making interactive (dynamic) graphs in R, also follows the grammar of graphics philosophy
library(DT)
library(scatterD3)

# Data -------
# for this part of the class you'll use your normalized and filtered data in log2 cpm
# make sure you have this object already in your work environment
# if you don't, go back to the Step2 script and generate it
head(log2.cpm.filtered.norm)
#give your data some more informative column names
colnames(log2.cpm.filtered.norm) <- sampleLabels
# Now we need to convert our datamatrix to a dataframe, while preserving the rownames as a new column in this dataframe
mydata.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneSymbol")
write_tsv(mydata.df, "normData.txt") # Note: this is the data you would use as input .gct file for GSEA analysis

# use dplyr 'mutate' function to add new columns based on existing data -------
mydata.df <- mutate(mydata.df,
                   uninfected.AVG = (uninf_rep1 + uninf_rep2 + uninf_rep1)/3, 
                   crypto.wt.AVG = (crypto.wt_rep1 + crypto.wt_rep2 + crypto.wt_rep3)/3,
                   crypto.mut.AVG = (crypto.mut_rep1 + crypto.mut_rep2 + crypto.mut_rep3)/3,
                   #now make columns comparing each of the averages above that you're interested in
                   LogFC.crypto.wt_vs_uninfected = (crypto.wt.AVG - uninfected.AVG),
                   LogFC.crypto.mut_vs_uninfected = (crypto.mut.AVG - uninfected.AVG))
#why is this type of approach to managing a spreadsheet useful?

#now look at this modified data table
mydata.df

# using dplyr 'arrange' and 'select' to sort your dataframe based on any variable -----
# first, we'll use dplyr "arrange" function to sort rows based on the values in a column of interest
# then we'll display 'select' only the columns we're interested in seeing
mydata.sort <- mydata.df %>%
  dplyr::arrange(desc(LogFC.crypto.wt_vs_uninfected)) %>% #note that this is the first time you've seen the 'pipe' operator
  dplyr::select(geneSymbol, LogFC.crypto.wt_vs_uninfected)


# use dplyr "filter" and "select" functions to pick out genes of interest  ----
#ways to tweek the 'select' function
#use : between two column names to select all columns between
#use 'contains', 'starts_with' or 'ends_with' to modify how you select
#can refer to columns using exact name or numerical indicator
#use boolean operators such as '&' (and), '|' (or), '==' (equal to), '!' (not)
mydata.filter <- mydata.df %>%
  dplyr::filter(geneSymbol=="MX1" | geneSymbol=="IRF1" | geneSymbol=="OAS2" | geneSymbol=="IRF3") %>%
  dplyr::select(geneSymbol, uninfected.AVG, crypto.wt.AVG)
mydata.filter

# you can also filter based on any regular expression
mydata.grep <- mydata.df %>%
  dplyr::filter(grepl('CXCL|IFI', geneSymbol)) %>%
  dplyr::select(geneSymbol, uninfected.AVG, crypto.wt.AVG)

mydata.grep

# make an interactive table ----
datatable(mydata.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'my cool table)',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(1:6), digits=3)

# make a simple scatter plot -----
ggplot(mydata.df, aes(x=uninfected.AVG, y=crypto.wt.AVG)) +
  geom_point(shape=3) +
  geom_point(size=1)

# make an interactive scatter plot ----
tooltip <- paste0("<b>","Symbol: ", mydata.df$geneSymbol, "</b><br>",
                  "<b>","untreated.AVG: ", "<b>", mydata.df$uninfected.AVG,
                  "<b>", "wt_crypto.AVG: ", "<b>", mydata.df$crypto.wt.AVG)

#plot the interactive graphic
scatterD3(mydata.df, x = uninfected.AVG, y = crypto.wt.AVG,
          lasso = TRUE,
          xlab = "uninfected.AVG", 
          ylab = "crypto.wt.AVG",
          colors = "red",
          point_opacity = 0.7,
          caption = list(title = "impact of infection"),
          tooltip_text = tooltip, hover_size = 3)

mydata.df %>% 
  ggvis(x= ~uninfected.AVG, y= ~crypto.wt.AVG, 
        key := ~geneSymbol,
        size := input_slider(10, 100),
        opacity := input_slider(0, 1)) %>% 
  add_tooltip(tooltip)
#notice that your R session is always active during an interactive session
