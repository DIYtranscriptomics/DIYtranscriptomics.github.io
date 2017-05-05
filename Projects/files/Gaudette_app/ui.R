library(shiny)
library(ggvis)
library(markdown)
shinyUI(fluidPage(
  headerPanel('Gene expression of plasma cell subsets'),
  sidebarPanel(
    selectInput("xvar", "X-axis variable", axis_vars, selected = "Splenic.Follicular_Bcell.B220P.AVG"),
    selectInput("yvar", "Y-axis variable", axis_vars, selected = "Bone_Marrow.BlimpP_PCs.B220N.AVG"),
    tags$small(paste0())
    ),
  fluidRow(
    includeMarkdown("include.md")),
  
  mainPanel(
    h2("Expression plot"),
    ggvisOutput("plot1")),
  
  fluidRow(
    column(6, DT::dataTableOutput('view'))
  )
))
