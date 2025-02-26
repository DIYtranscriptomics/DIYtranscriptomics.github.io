# Introduction ----
# The goal of this script is to analyze data from LEMIS, tracking the import of animals and animal products into US ports.
# We'll use a combination of tidyverse and ggplot2 to explore the data and answer a series of questions.
#
# Load tidyverse and plotly packages ----


# read in lemis data, stored locally as lemis_cleaned.tsv ----


# let's experiment with the prompt above.  What happens if we change .tsv to .csv...or remove these details altogether?

# take a look at the unique levels for any of the variables in this dataset


# Question 1: Use the filter function from dplyr to find the most common (by 'quantity') live mammal taken from the wild for import into the US. ----
lemis.live <- lemis %>%
  dplyr::filter(description == "Live specimens (live animals or plants)") %>%
  dplyr::filter(class == "Mammalia") %>%
  dplyr::filter(source == "Specimens taken from the wild") %>%
  dplyr::arrange(desc(quantity))

# Question 2: Building on your analysis above, filter the dataset to get live animals (use 'generic_name') imported for the purposes of "scientific" or "Biomedical research"


# plot this lemis.live data you generated above using ggplot2


# generate the plot above again, but this time using ploty


# Question 3: From which countries do we import macaques? ----


# Question 4: From which countries from which we import live bats? ----


# Question 5: For what 'purpose' are these imported bats used? ----


# Question 6: How does the type of bat (use 'specific_name') imported differ between countries (use 'facet_wrap' in your ggplot code)? ----


# Question 7: Identify the most expensive (by ‘value’) shipment of "Live specimens (live animals or plants)" to enter the US? ----


# Question 8: Identify the most expensive shipment of any kind (live or not)? ----



