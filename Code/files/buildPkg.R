# Introduction -----
# this script walks through the basic steps involved in turning a project into an R package

# Load packages ----
library(devtools) #lots of tools for developers
library(roxygen2) #for easily creating package documentation 
library(usethis) #automates many steps involved in creating and managing packages and projects in R
library(available) #let's you check if a package name is available. Download using: devtools::install_github("ropenscilabs/available")

# Construct a basic package ----
# you can check to see if you package name is available, using the available function
available("DIYtranscriptomics", browse = FALSE)

# running the 'create_package()' function will create an R package and open a new Rsession in that package
create_package("~/Desktop/DIYtranscriptomics") 

# Add fuction scripts to 'R' directory----
# functions are the heart of any package
# To begin, it's a good idea to use a single R script for each function
# first, change your working directory to be the new R package/project you just set-up
setwd("~/Desktop/DIYtranscriptomics")
# the 'use_r()' function will create a new R script in /R with the name you provide
# the empty script will open in RStudio and you can paste your function code in and save it
use_r("DIYprofile")

# Add data to your package ----
use_data(targets) #will create a data directory with your targets saved as .rdata file
use_data(myDGEList) #adds myDGEList.rdata to the data directory

# Document scripts and objects ----
# now that you've added data and functions, 
# you need to document them so that people who download your package will understand what these are and how they relate to your analysis project
# It is important to add documentation to your data and functions
# we do this using the roxygen package
document() #creates an empty 'man' folder in your package directory
# add comments to the top of any function.R script using "#'" comment tags
# see examples at: https://r-pkgs.org/man.html
# then run document() again
document()
# now let's take a look at how this looks for own R function using '?' to view our own help docs!
?DIYprofile
# we should also document our data, but this happens a bit differently
# we can't really add roxygen comments to the top of an R data object
# so instead we add a 'data.R' script to /data and add roxygen comments there for each dataset
use_r("data")
document()
# take a look at what this created
?myDGEList
?targets

# Add readme and vignette ----
# running 'use_readme_md()' will create a markdown readme and open it in a new tab of Rstudio
# edit this readme to describe your package. It will be the first thing people see when they look at your package on GitHub
use_readme_md()
# 'use_vignette()' will create a vignette for your package and open in a new tab
# copy and past your Rmarkdown content into this vignette
use_vignette(name = "DIYtranscriptomics") 

# Set dependencies for your package ----
# This is where you decide which packages will be downloaded/installed when someone installs your package
# see here for more info: https://r-pkgs.org/description.html#dependencies
# running 'use_package()' will add a package to the imports field of your DESCRIPTION file
use_package("ggplot2")
use_package("tibble")
use_package("tidyr")
use_package("edgeR")
use_package("rmarkdown")
use_package("tinytex")
use_package("knitr")
use_package("dplyr")
use_package("Biobase")
use_package("DT")
use_package("GSEABase")
use_package("GSVA")
use_package("RColorBrewer")
use_package("clusterProfiler")
use_package("cowplot")
use_package("enrichplot")
use_package("gplots")
use_package("gprofiler2")
use_package("gt")
use_package("limma")
use_package("matrixStats")
use_package("msigdbr")
use_package("plotly")

# Check and edit your package ----
# switch to your RStudio project for the package you've created
# take a look at your DESCRIPTION file. You should see all of the packages above in the 'Imports' field
# edit the 'Description' field of the DESCRIPTION file to suit you
# edit the 'Authors' field as well

# Now give your package a license
use_mit_license("Daniel Beiting")

# We're ready to check the package to see if it works!
# can run the 'check()' function below, or use the 'build' menu in RStudio
# either way, running 'check' will produce an output folder with your Rmarkdown in the 'docs' directory
# load all R code from /R using 'load_all()' and test that your function works
load_all() # probably best to run in the RStudio session for your new package
check() # probably best to run in the RStudio session for your new package
# If check proceeded with no errors or warning, then you are ready to build your package!
build() # probably best to run in the RStudio session for your new package

# Autenticate to use GitHub ----
# You only need to do this authentication once
use_git_config(user.name = "dpbisme", user.email = "danielbeiting@gmail.com")
git_sitrep()
# now get a personal access token (PAT) so you can create a github repo directly from RStudio
# running 'browse_github_token' will open your browser and navigate to the github PAT page
# scroll to the bottom of this page and and click 'get token'
# be SURE to copy this token. Once you close the webpage you will never see it again.
browse_github_token() 
# run 'edit_r_environment' and paste in your PAT as a new line in the file that .renviron file that opens
# line should read 'GITHUB_PAT=...', where ... is your PAT
edit_r_environ() #you'll need to restart R after this by going to RStudio Session menu -> restart R
# since you restarted R, you'll need to reload devtools and roxygen from the top of this script

# Put your package on GitHub ----
use_git() #makes your project version controlled
use_github() #pushes contents up to github

# Anyone can reproduce your analysis by: ----
# 1. downloading your package repo from github
# 2. opening the .rproj file in the package to launch RStudio
# 3. running devtools::build() in RStudio
# 4. Opening the resulting .tar.gz package that is created by build()
# 5. finding the Rmarkdown html file in the /inst/docs directory in the unzipped project
# 6. If so desired, they could also run devtools::install() to install and load your library in their environment (which gives them access to your functions)
