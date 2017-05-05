#this is the script that sets up the shiny server with code necessary to generate the correct graphs, data and reactives
#this script begins by loading an R object that I've called 'myData'
#in my case the 'myData' object was saved from a previous R analysis and is simply a customized dataframe of gene expression data produced using dplyr
#for an example of a 'myTPM.df' object, see my 'Step3_dataWrangling.R' script where I use the mutate command from dplyr
#to customize the script for your own data (could be RNAseq, array, microbiome, ...really anything), you'll need to do the following:
#1. load your own data on line 14 below, instead of 'myTPM.df'
#2. change the 'tooltip' (line 17) to be something you'd like to see when you mouse over data points on the final plot (in my case, it's set to show gene symbol)
#3. I've used the 'ggvis' package below to build a scatter plot for looking at pairwise gene expression.  You may want to change this as well
#4. Also replace the 'myTPM.df' piped to the ggvis function (line 24) and then called by the 'renderDataTable' function (line 35) at the end so that you build a reactive datatable from your own data

library(shiny)
library(DT)

#using 'myData' object created from dplyr mutate command
load("myTPM.df")

shinyServer(function(input, output) {
  tooltip <- function(data, ...) {
    paste0("<b>", data$geneSymbol, "</b>")
  }
  
  vis <- reactive({
    # reactive labels for axes of the scatterplot
    xvar <- prop("x", as.symbol(input$xvar))
    yvar <- prop("y", as.symbol(input$yvar))
   #passing 'myTPM.df' to ggvis to make a graph 
    myTPM.df %>% 
      ggvis(x= xvar, y= yvar, key := ~geneSymbol) %>% 
      add_tooltip(tooltip)
  }) 
  
  vis %>% bind_shiny("plot1")


  output$view <- DT::renderDataTable(myTPM.df, server = TRUE)
})

