install.packages("openxlsx", dependencies=TRUE)
install.packages("ggplot2", dependencies=TRUE)
install.packages("plotly", dependencies=TRUE)
install.packages("treemapify", dependencies=TRUE)

#########################brighter colours*****subtle colours(below)########################

library(tidyverse)
library(ggfittext)
library(ggplot2) 
library(treemapify)

treemap_plot <- read.csv("~/Documents/PhD/Chapter_2/amr_files/perl_flago_cov/Overall_flago.csv")

ggplot(treemap_plot, aes(area = value, fill = parent, colour=id, label=id, subgroup=parent)) +
  geom_treemap(fill = "black") +
  geom_treemap(aes(alpha = value)) +
  geom_treemap_subgroup_border(colour="white") +
  geom_treemap_text(fontface = "italic",
                    colour = "white",
                    place = "centre",
                    grow = F,
                    reflow=T) +
  geom_treemap_subgroup_text(place = "centre",
                             grow = T,
                             alpha = 0.5,
                             colour = "#FAFAFA",
                             min.size = 0) +
  scale_alpha_continuous(range = c(0.4, 1))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))


#########################subtle colours########################
library(dplyr)
library(colorspace)  # install via: install.packages("colorspace", repos = "http://R-Forge.R-project.org")
library(scales)

# number of palettes needed
n <- length(unique(treemap_plot$parent))

# now calculate the colors for each data point
df2 <- treemap_plot %>%
  mutate(index = as.numeric(factor(parent))- 1) %>%
  group_by(index) %>%
  mutate(
    max_value = max(value),
    color = gradient_n_pal(
      sequential_hcl(
        6,
        h = 360 * index[1]/n,
        c = c(45, 20),
        l = c(30, 80),
        power = .5)
    )(value/max_value)
  )

ggplot(df2, aes(area = value, fill = color, label=id, subgroup=parent)) +
  geom_treemap() +
  geom_treemap_subgroup_border(colour="white") +
  geom_treemap_text(fontface = "italic",
                    colour = "white",
                    place = "centre",
                    grow = F,
                    reflow=T) +
  geom_treemap_subgroup_text(place = "centre",
                             grow = T,
                             alpha = 0.5,
                             colour = "#FAFAFA",
                             min.size = 0) +
  scale_fill_identity()