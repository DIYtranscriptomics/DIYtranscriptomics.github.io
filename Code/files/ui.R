library(shiny)
library(ggvis)
library(markdown)
shinyUI(fluidPage(theme="cerulean.css",
  headerPanel('Host cell response to infection with Cryptosporidium'),
  #text that describes the experiment, contained in a separate file called 'include.md'
  #headerPanel(
    #img(src='PennVet_Logo.jpg', align = "right")),
  headerPanel(
    includeMarkdown("Abstract.md")),
  #set the dropdown menus used to select axes
  sidebarPanel(
    selectInput("xvar", "X-axis variable", axis_vars, selected = "untreated.AVG"),
    selectInput("yvar", "Y-axis variable", axis_vars, selected = "wt_crypto.AVG"),
    tags$small(paste0()),
    includeMarkdown("Instructions.md")
    ),
  #now the graph
  sidebarPanel(width=8,
    h2("Expression plot"),
    ggvisOutput("plot1")),
  
  #table of data
  fluidRow(
    column(6, DT::dataTableOutput('view'))
  )
))
