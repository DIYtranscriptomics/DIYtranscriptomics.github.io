library(shiny)
library(ggvis)
library(markdown)
shinyUI(fluidPage(
  headerPanel('Role of Caspase-8 on TLR Signaling'),
  #text that describes the experiment, contained in a separate file called 'include.md'
  #headerPanel(
    #img(src='PennVet_Logo.jpg', align = "right")),
  headerPanel(
    includeMarkdown("Abstract.Rmd")),
  #set the dropdown menus used to select axes
  sidebarPanel(
    selectInput("xvar", "X-axis variable", axis_vars, selected = "B6.untreated.1"),
    selectInput("yvar", "Y-axis variable", axis_vars, selected = "Ripk3.untreated.1"),
    tags$small(paste0()),
    includeMarkdown("Instructions.Rmd")
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
