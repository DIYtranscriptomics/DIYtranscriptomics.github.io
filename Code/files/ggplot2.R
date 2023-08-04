# load packages and data ----
library(tidyverse) # you're already familiar with this from module 3
library(RColorBrewer) # a great R package for working with colors
library(plotly) # makes it easy to convert static plots to interactive
# recall from module 1 that we have lots of datasets prepackaged with R
# let's take a look at this list and scroll down to the tidyverse datasets
data()
# take a look at the mpg data
mpg

# let's iteratively build a plot with ggplot2
ggplot() # nothing
ggplot(mpg) # even if we provide some data, we still don't get a plot

ggplot(mpg) + # the '+' is an important part of ggplot b/c it allows layering!
  aes(x=hwy, y=cty) # with data + mappings...still nothing.  What are we missing?

ggplot(mpg) +
  aes(x=hwy, y=cty) +
  geom_point() # what is a 'geom'?  

ggplot(mpg) + 
  aes(x=hwy, y=cty) +
  geom_point(aes(size = cyl)) # we can map more of the data to point color

ggplot(mpg) +
  aes(x=hwy, y=cty) +
  geom_point(aes(size = cyl)) +
  geom_smooth(method ="lm") # additional plots can be overlayed

ggplot(mpg) +
  aes(x=hwy, y=cty) +
  geom_point(aes(size = cyl)) +
  geom_smooth(method ="lm") +
  # let's add some more human readable axis labels and plot titles
  labs(x="mpg on highway", y="mpg in city",
       title = "Analysis of MPG dataset",
       subtitle = "examining gas efficiency by year and engine type")

ggplot(mpg) +
  aes(x=hwy, y=cty) +
  geom_point(aes(size = cyl)) +
  geom_smooth(method ="lm") +
  labs(x="mpg on highway", y="mpg in city",
       title = "Analysis of MPG dataset",
       subtitle = "examining gas efficiency by year and engine type") +
  theme_classic() # themes make everything look more polished

# the ggsave function saves the last plot produced to a file
ggsave("plot.png", width = 8, height = 5)

# it is easy to facet plots based on any available data
ggplot(mpg) +
  aes(x=hwy, y=cty) +
  geom_point(aes(size = cyl)) +
  geom_smooth(method ="lm") +
  labs(x="mpg on highway", y="mpg in city",
       title = "Analysis of MPG dataset",
       subtitle = "examining gas efficiency by year and engine type") +
  theme_bw() +
  facet_wrap(~year) # this is all that is needed to facet

# ggplots can be stored as an object
myPlot <- ggplot(mpg) +
  aes(x=hwy, y=cty) +
  geom_point(aes(size = cyl)) +
  geom_smooth(method ="lm") +
  labs(x="mpg on highway", y="mpg in city",
       title = "Analysis of MPG dataset",
       subtitle = "examining gas efficiency by year and engine type") +
  theme_bw() +
  facet_wrap(~year) # this is all that is needed to facet

# only takes one line of code to make this plot interactive using plotly
ggplotly(myPlot)

# working with colors in R
#Some useful examples: rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]).
display.brewer.all() # note that there are qualitative, diverging, and sequential palletes to choose from
display.brewer.all(colorblindFriendly = TRUE)
brewer.pal(5,'RdYlBu') # notice that this produces a simple character vector of hexcodes
# RStudio recognizes hexcodes and produces the color directly in the browser

# Recreating Minard's 'Napoleon's March' in ggplot2 ----
# code from Andrew Heiss: https://www.andrewheiss.com/blog/2017/08/10/exploring-minards-1812-plot-with-ggplot2/
library(tidyverse)
library(ggmap)
library(ggrepel)
library(gridExtra)
library(pander)

# reading in data (came from here: https://github.com/andrewheiss/fancy-minard/tree/master/input/minard)
cities <- read.table("cities.txt",
                     header = TRUE, stringsAsFactors = FALSE)

troops <- read.table("troops.txt",
                     header = TRUE, stringsAsFactors = FALSE)

temps <- read.table("temps.txt",
                    header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(date = dmy(date))  # Convert string to actual date

temps.nice <- temps %>%
  mutate(nice.label = paste0(temp, "°, ", month, ". ", day))

march.1812.plot.simple <- ggplot() +
  geom_path(data = troops, aes(x = long, y = lat, group = group,
                               color = direction, linewidth = survivors),
            lineend = "round") +
  geom_point(data = cities, aes(x = long, y = lat),
             color = "#DC5B44") +
  geom_text_repel(data = cities, aes(x = long, y = lat, label = city),
                  color = "#DC5B44", family = "Open Sans Condensed Bold") +
  scale_size(range = c(0.5, 10)) +
  scale_colour_manual(values = c("#DFC17E", "#252523")) +
  guides(color = FALSE, size = FALSE) +
  theme_nothing()

# Change the x-axis limits to match the simple map
temps.1812.plot <- ggplot(data = temps.nice, aes(x = long, y = temp)) +
  geom_line() +
  geom_label(aes(label = nice.label),
             family = "Open Sans Condensed Bold", size = 2.5) +
  labs(x = NULL, y = "° Celsius") +
  scale_x_continuous(limits = ggplot_build(march.1812.plot.simple)$layout$panel_ranges[[1]]$x.range) +
  scale_y_continuous(position = "right") +
  coord_cartesian(ylim = c(-35, 5)) +  # Add some space above/below
  theme_bw(base_family = "Open Sans Condensed Light") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks = element_blank(),
        panel.border = element_blank())

# Combine the two plots
both.1812.plot.simple <- rbind(ggplotGrob(march.1812.plot.simple),
                               ggplotGrob(temps.1812.plot))

# Adjust panels
panels <- both.1812.plot.simple$layout$t[grep("panel", both.1812.plot.simple$layout$name)]

# Because this plot doesn't use coord_equal, since it's not a map, we can use
# whatever relative numbers we want, like a 3:1 ratio
both.1812.plot.simple$heights[panels] <- unit(c(3, 1), "null")

grid::grid.newpage()
grid::grid.draw(both.1812.plot.simple)


# Recreating the gapminder visualization in ggplot2 ----
library(ggplot2)
library(gganimate)
library(gapminder)

ggplot(gapminder) +
  aes(x=gdpPercap, y=lifeExp, size = pop, color = continent, frame = year) +
  geom_point() +
  scale_x_log10() +
  theme_bw() +
  transition_time(year) + # this line and the next are possible because of gganimate
  ease_aes('linear')

# Recreating Ed Hawkins' 'Climate Spiral' in ggplot2 ----
# code from Pat Schloss: https://github.com/riffomonas/climate_viz
library(tidyverse)
library(gganimate)

t_diff <- read_csv("GLB.Ts+dSST.csv", skip = 1, na = "***") %>%
  select(year = Year, month.abb) %>%
  pivot_longer(-year, names_to="month", values_to="t_diff") %>%
  drop_na()

# last_dec <- t_diff %>%
#   filter(month == "Dec") %>%
#   mutate(year = year + 1,
#          month = "last_Dec")

next_jan <- t_diff %>%
  filter(month == "Jan") %>%
  mutate(year = year - 1,
         month = "next_Jan")

t_data <- bind_rows(t_diff, next_jan) %>%
  mutate(month = factor(month, levels = c(month.abb, "next_Jan")),
         month_number = as.numeric(month)) %>%
  arrange(year, month) %>%
  filter(year != 1879) %>%
  mutate(step_number = 1:nrow(.))

annotation <- t_data %>%
  slice_max(year) %>%
  slice_max(month_number)

temp_lines <- tibble(x = 12, y = c(1.5, 2.0), labels = c("1.5\u00B0C", "2.0\u00B0C"))
month_labels <- tibble(x = 1:12, labels = month.abb, y = 2.7)

a <- t_data %>% 
  ggplot(aes(x=month_number, y=t_diff, group=year, color=year)) +
  geom_rect(aes(xmin=1, xmax=13, ymin=-2, ymax=2.4),
            color="black", fill="black",
            inherit.aes = FALSE) +
  geom_hline(yintercept = c(1.5, 2.0), color="red") +
  geom_label(data = temp_lines, aes(x=x, y=y, label=labels),
             color = "red", fill = "black", label.size = 0,
             inherit.aes=FALSE) +
  geom_text(data = month_labels, aes(x=x, y=y, label = labels),
            inherit.aes = FALSE, color="white",
            angle = seq(360 - 360/12, 0, length.out = 12)) +
  geom_label(aes(x = 1, y=-1.3, label = year),
             color="white", fill="black",
             label.padding = unit(50, "pt"), label.size = 0,
             size=6) +
  geom_line() +
  scale_x_continuous(breaks=1:12,
                     labels=month.abb, expand = c(0,0),
                     sec.axis = dup_axis(name = NULL, labels=NULL)) +
  scale_y_continuous(breaks = seq(-2, 2, 0.2),
                     limits = c(-2, 2.7), expand = c(0, -0.7), 
                     sec.axis = dup_axis(name = NULL, labels=NULL)) +
  scale_color_viridis_c(breaks = seq(1880, 2020, 20),
                        guide = "none") +
  coord_polar(start = 2*pi/12) +
  labs(x = NULL,
       y = NULL,
       title = "Global temperature change (1880-2022)") +
  theme(
    panel.background = element_rect(fill="#444444", linewidth=1),
    plot.background = element_rect(fill = "#444444", color="#444444"),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(color="white", size=13),
    plot.title = element_text(color="white", hjust = 0.5,size = 15)) +
  transition_manual(frames = year, cumulative = TRUE)

animate(a, width=4.155, height=4.5, unit="in", res=300,
        renderer = av_renderer("climate_spiral.mp4"))


