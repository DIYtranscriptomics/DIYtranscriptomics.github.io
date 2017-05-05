my.df <- read.delim("data.unfiltered.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE)
head(my.df)
dim(my.df)
class(my.df)

#load dplyr for doing most of the heavy lifting in this hackdash
library(dplyr)

#I know there are 4 main 'verbs' in dplyr
#1. select - picks columns
#2. arrange - sort rows
#3. filter - filter rows
#4. mutate - create new columns based on existing data in the spreadsheet

# begin by selecting all columns with names that contain a 'F' (will grab all samples from female worms)
my.df.subset  <- dplyr::select(my.df, contains("F"), description, transcript_id)
head(my.df.subset)

# now use the '-contains' option in dplyr to drop all the samples from LEPZQ and NMRI strain parasites
# this leaves us with just Female LE strain worms
my.df.subset <- dplyr::select(my.df.subset, -contains("LEPZQ"))
my.df.subset <- dplyr::select(my.df.subset, -contains("NMRI"))
head(my.df.subset)

# now use the 'mutate' function to create new AVG and LogFC columns
my.df.subset.mod <- mutate(my.df.subset,
                           control.AVG = (FCtl_LE_1 + FCtl_LE_2 + FCtl_LE_3)/3,
                           tp3hr.AVG = (F3h_LE_1 + F3h_LE_2 + F3h_LE_3)/3,
                           tp12hr.AVG = (F12h_LE_1 + F12h_LE_2 + F12h_LE_3)/3,
                           tp24hr.AVG = (F24h_LE_1 + F24h_LE_2 + F24h_LE_3)/3,
                           LogFC_3hr.vs.ctrl = tp3hr.AVG - control.AVG,
                           LogFC_12hr.vs.ctrl = tp12hr.AVG - control.AVG,
                           LogFC_24hr.vs.ctrl = tp24hr.AVG - control.AVG,
                           transcript_id = transcript_id,
                           description = description)
                           
head(my.df.subset.mod)

# arrange this new dataframe in descending order based on the LogFC column for 24hr vs control
my.df.subset.mod <- arrange(my.df.subset.mod, desc(LogFC_24hr.vs.ctrl))
head(my.df.subset.mod)
dim(my.df.subset.mod)

# filter the rows to leave only genes with 2fold or greater induction at 24hr
my.df.subset.FC.select <- filter(my.df.subset.mod, LogFC_24hr.vs.ctrl >= 2)
dim(my.df.subset.FC.select)

# using this filtered dataframe, select just the AVG, LogFC columns, and annotation columns
my.df.subset.FC.select2 <- select(my.df.subset.FC.select, contains("AVG"), contains("LogFC"), description, transcript_id)

head(my.df.subset.FC.select2)
dim(my.df.subset.FC.select2)

# end by creating an interactive table with the result from above
library(DT)
datatable(my.df.subset.FC.select2, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'I win the hackdash!',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(1:6), digits=3)
