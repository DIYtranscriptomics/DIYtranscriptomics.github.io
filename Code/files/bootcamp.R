# Introduction----
# good R scripts should be extensively annotated with comments like this one! 
# comments should explain what a particular piece of code is trying to accomplish and why
# commenting with '#' is also a handy way to toggle on or off specific lines of code
# Before getting started, be sure to set your working directory and make sure all files you need are in that directory

# Install and load packages----
# install tidyverse
install.packages("tidyverse") 
# load tidyverse
library(tidyverse)
# note that R 'looks' for packages to download from specific repositories
setRepositories() # note that this is the first time we've used a function with no 'arguments'

# Getting data into and out of R----
# R comes prepackaged with many example datasets, which you can browse with:
data()
# and then view any of these internal datasets with:
View(Titanic)
# a handy way to view data is to use the glimpse function from dplyr
glimpse(Titanic)
# read a data file in from your computer using the read_delim function from the readr package 
my.dataframe <- read_delim("swine_study.txt") # first time you've seen the assignment operator
# write a R object out to a file
write_tsv(my.dataframe,"test.txt")

# Exploring common data structures in R ----
# vector
myNumericVector <- c(2,4,8,3,1) #notice I don't use spaces in variable names!
my_character_vector <- c("have", "a", "very", "very", "good", "day!")
# one great thing about R is that it allows 'vectorization', where functions act on all elements of a data structure.
myNumericVector + 10
# another handy thing we can do with vectors is to factor them.  This allows R to see 'levels' for each element in the vector
my_character_vector <- factor(my_character_vector)
# matrix
my.matrix <- as.matrix(read_delim("matrix.txt")) # notice how we can nest functions
# dataframe
my.dataframe # we already created this object in the section above
# we can access any column of our dataframe using the '$' operator
my.dataframe$SampleID
# list
myList <- c(my.dataframe, myNumericVector, my_character_vector)
# we can evaluate data structures with logical expressions
myLogical <- is.numeric(myNumericVector)
myLogical <- is.numeric(my_character_vector)
# you can always check what type of data an object is either by looking in the Environment browser, or by:
class(myNumericVector)

# Subsetting data structures in R----
myVectorSubset <- my_character_vector[1:3]
myMatrixSubset <- my.matrix[,]
myDataframeSubset <- my.dataframe[,1:4]
myListSubset <- myList[1]
# we also use expressions to evaluate and subset our dataframe
my.dataframe$current_antibiotics == "Yes"
no_antibiotics <- my.dataframe[my.dataframe$current_antibiotics == "Yes", ]

# Finding help directly within R----
?read_delim # typing '?' followed by a function name will display help docs for that function
??delim # typing '??' followed by a keyword will search ALL help docs for that word
# note that you can highlight and run R code directly from the help docs!

# Other useful stuff----
# it's easy to turn any comment into a code chunk by adding '----' and end of line
