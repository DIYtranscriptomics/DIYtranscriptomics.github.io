# PREPARE: read data into R and explore ----
library(tidyverse)
lemis <- read_delim("lemis_cleaned.tsv")
#take a look at all the variables and values ('levels' of a variable) in this dataset
glimpse(lemis)
#here's how you can see all the levels of any variable in this dataset
unique(lemis.live$description)

# QUESTION_01: Identify the most common (by 'quantity') live mammal taken from the wild for import into the US. ----
# note that I will continue to use this 'live.lemis' dataframe many times throughout this script
lemis.live <- lemis %>%
  dplyr::filter(description == "Live specimens (live animals or plants)") %>%
  dplyr::filter(class == "Mammalia") %>%
  dplyr::filter(source == "Specimens taken from the wild") %>%
  dplyr::arrange(desc(quantity))

# note that the above solution ranks by the largest shipment, but a smarter way to do with would consider the sum of all shipments for each animal type.
# this is pretty easy to do using 'group_by' and 'summarize' from dplyr
lemis.live.summary <- lemis %>%
  dplyr::filter(description == "Live specimens (live animals or plants)") %>%
  dplyr::filter(class == "Mammalia") %>%
  dplyr::filter(source == "Specimens taken from the wild") %>%
  group_by(generic_name) %>%
  summarize(sum = sum(quantity, na.rm = TRUE)) %>%
  dplyr::arrange(desc(sum))
# bears still win :)


# QUESTION_02: Building on your analysis above, produce a plot showing live animals (use 'generic_name') imported for the purposes of science/research ----
lemis.science <- lemis %>%
  dplyr::filter(description == "Live specimens (live animals or plants)") %>%
  dplyr::filter(class == "Mammalia") %>%
  dplyr::filter(source == "Specimens taken from the wild") %>%
  dplyr::filter(purpose == "Scientific" | purpose == "Biomedical research")

ggplot(lemis.science) +
  aes(y=generic_name, x=quantity) +
  geom_col() +
  theme_bw()

# QUESTION_03: Identify the countries from which we import the most macaques (again, use a simple plot). ----

lemis.live.macaques <- lemis.live %>%
  dplyr::filter(generic_name == "MACAQUE")

ggplot(lemis.live.macaques) +
  aes(y=country_origin, x=quantity) +
  geom_col()
#how does this result change if we don't require them to come from the wild?

# QUESTION_04 and _05: Using the same approach as above, create a plot showing the countries from which we import live bats. ----
lemis.live.bats <- lemis.live %>%
  dplyr::filter(generic_name == "BAT")

ggplot(lemis.live.bats) +
  aes(y=country_origin, x=quantity) + #see why these bats are imported by changing to 'y=purpose'
  geom_col()

# QUESTION_06: How does the type of bat (use 'specific_name') imported differ between countries (hint: use ```facet_wrap``` in your ggplot code)? ----
ggplot(lemis.live.bats) +
  aes(y=specific_name, x=quantity) +
  geom_col() +
  facet_wrap(~country_origin)

# QUESTION_07: Identify the most expensive (by 'value') shipment of live mammals to enter the US ----
lemis.live.expensive <- lemis.live %>%
  dplyr::top_n(1, value) # note the use of top_n here

# QUESTION_08: How does the answer above compare with the most expensive shipment of any kind ----
lemis.all.expensive <- lemis %>%
  dplyr::top_n(1, value)

# QUESTION_09: You are alerted to a concerning new viral disease of humans that is believed to originate from Fruit bats.  Identify the US city that would would be most likely to be exposed to such a virus from the import of live fruit bats. ----
lemis.bats.fruit <- lemis.live %>%
  dplyr::filter(specific_name == "FRUIT" | specific_name == "TALAUD FRUIT" | specific_name == "RED-NECKED FRUIT")

# alternative approach using grepl
lemis.bats.fruit <- lemis.live %>%
  dplyr::filter(grepl('FRUIT', specific_name))

ggplot(lemis.bats.fruit) +
  aes(y=port, x=quantity) +
  geom_col()

# QUESTION_10: A recent case of Anthrax in NYC was traced back to a contaminated Wildebeest hide that was stretched and used to make a traditional drum.  Through which port(s) did this animal product most likely enter the country? ----
lemis.skin <- lemis %>%
  dplyr::filter(description == "Skin piece (raw or tanned including scraps)") %>%
  dplyr::filter(generic_name == "WILDEBEEST")

ggplot(lemis.skin) +
  aes(y=port, x=quantity) +
  geom_col()




# BONUS IDEAS ----
# Generate a report on confiscated goods, here's some code to get you started
lemis.seized <- lemis %>%
  dplyr::filter(disposition == "Seized") %>%
  dplyr::arrange(desc(quantity)) %>%
  dplyr::top_n(25,quantity)

ggplot(lemis.seized) +
  aes(y=generic_name, x=quantity) +
  geom_col()
# Use ggAnimate to create a animation of import data over time
library(gganimate)
library(gifski) #only need to download once and gganimate will load
library(av) #only need to download once and gganimate will load
# restart R after installing the gifski and av packages above, then you're ready to go
# let's use the data for seized shipments that we generated above
ggplot(lemis.seized) +
  aes(y=port, x=quantity, fill=generic_name) +
  geom_col(position = "identity") +
  # Here comes the gganimate specific bits
  labs(title = 'Year: {frame_time}') +
  transition_time(as.integer(shipment_year)) +
  ease_aes('linear')


# to save a gif produced by gganimate above, store your ggplot object (e.g., as 'plot')
# then pass to 'anim_save("filename.gif", plot)'




