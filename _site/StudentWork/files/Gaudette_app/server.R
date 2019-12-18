#this is the script that sets up the shiny server with code necessary to generate the correct graphs, data and reactives
#this script begins by loading an R object that I've called 'myData'
#in my case the 'myData' object was saved from a previous R analysis and is simply a customized dataframe of gene expression data produced using dplyr
#for an example of a 'myData' object, see my 'Step3_dataExploration_part2.R' script where I use the mutate command from dplyr
#to customize the script for your own data (could be RNAseq, array, microbiome, ...really anything), you'll need to do the following:
#1. load your own data instead of 'myData'
#2. change the 'tooltip' to be something you'd like to see when you mouse over data points on the end plot (in my case, it's set to show gene symbol)
#3. I've used the 'ggvis' package below to build a scatter plot for looking at pairwise gene expression.  You may want to change this as well
#4. Also replace the 'myData' called by the 'renderDataTable' function at the end so that you build a reactive datatable from your own data

library(shiny)
library(DT)

load("myData")
myDataAVG <- myData[,-1:-18]

shinyServer(function(input, output, session) {
  tooltip <- function(data, ...) {
    paste0("<b>","Symbol: ", data$geneSymbols, "</b>")
  }
  
  vis <- reactive({
    # Lables for axes
    xvar_name <- names(axis_vars)[axis_vars == input$xvar]
    yvar_name <- names(axis_vars)[axis_vars == input$yvar]
    xvar <- prop("x", as.symbol(input$xvar))
    yvar <- prop("y", as.symbol(input$yvar))
    myData %>% 
      ggvis(x= xvar, y= yvar, key := ~geneSymbols) %>% 
      add_tooltip(tooltip)
  }) 
  
  vis %>% bind_shiny("plot1")


  output$view <- DT::renderDataTable(myDataAVG, server = TRUE)
})

