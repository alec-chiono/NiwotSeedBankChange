# Vizualize raw counts as bar plot

# PACKAGES ----
library(tidyverse)
library(ggbreak)
library(ggthemes)
library(ggimage)
library(patchwork)

# DATA ----
# Source data download function
source("scripts/00_DownloadData.R")

# Download and wrangle data
seed_count_df <-
  download_data(seed = TRUE)$seedbank_composition.ac_hh.data %>%
  filter(
    substr(USDA_code, 1, 1) != 2 & USDA_code != "CAREX" & USDA_code != "POA"
  ) %>% #remove records not identified to species
  group_by(USDA_name, year) %>%
  summarize(sum = sum(count), .groups = "drop") %>%
  mutate(
    Species = sub("\\s*var\\..*$", "", USDA_name),
    year = factor(year)
  )

# VIZ ----
# Growth form icon files
gram <- "images/graminoid.png"
cush <- "images/cushion.png"
forb <- "images/forb.png"

# Assign species growth form
symbol_map <- c(
  "Allium geyeri" = forb,
  "Androsace septentrionalis" = forb,
  "Arenaria fendleri" = forb,
  "Calamagrostis purpurascens" = gram,
  "Campanula rotundifolia" = forb,
  "Deschampsia cespitosa" = gram,
  "Draba aurea" = forb,
  "Draba streptocarpa" = forb,
  "Erysimum capitatum" = forb,
  "Festuca brachyphylla" = gram,
  "Festuca minutiflora" = gram,
  "Geum rossii" = forb,
  "Juncus drummondii" = gram,
  "Minuartia obtusiloba" = cush,
  "Minuartia rubella" = forb,
  "Noccaea montana" = forb,
  "Phlox pulvinata" = cush,
  "Poa alpina" = gram,
  "Polygonum viviparum" = forb,
  "Saxifraga rhomboidea" = forb,
  "Sedum lanceolatum" = forb,
  "Silene acaulis" = cush
)

# Order Species names by highest->lowest count in 1989
species_order <-
  seed_count_df %>%
  filter(year == 1989) %>%
  arrange(sum, desc = TRUE) %>%
  pull(Species)

# Data.frame to plot growth form icons under species names
icon_df <-
  tibble(
    Species = names(symbol_map),
    image = unname(symbol_map),
    y_icon = -2.2
  ) %>%
  mutate(Species = factor(Species, levels = species_order))

# Main bar plot
p_main <- seed_count_df %>%
  mutate(Species = factor(Species, levels = species_order)) %>%
  ggplot(aes(x = Species, y = sum, fill = year, colour = year)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = -0.158, color = "black", linewidth = 0.4) +
  geom_segment(
    x = rep(names(symbol_map), 2),
    y = -0.16,
    yend = -0.75,
    color = "black",
    linewidth = 0.4
  ) +
  geom_image(
    data = icon_df,
    aes(x = Species, y = y_icon, image = image),
    inherit.aes = FALSE,
    size = 0.065,
    by = "width"
  ) +
  scale_fill_manual(
    values = c("1989" = "skyblue3", "2023" = "red4"),
    guide = "none"
  ) +
  scale_color_manual(
    values = c("1989" = "skyblue3", "2023" = "red4"),
    guide = "none"
  ) +
  scale_y_continuous(
    name = "# Germinable Seeds",
    expand = expansion(mult = c(0, 0.05), add = c(3.2, 0)),
    #breaks = c(0, 5, 10, 15, 20, 80, 90, 100),
    limits = c(NA, 100)
  ) +
  #scale_x_discrete(expand = expansion(add = c(2.5, 0.2))) +
  scale_y_break(c(21, 80), scales = 0.25) +
  coord_cartesian(clip = "off") +
  theme_tufte() +
  theme(
    axis.text.x = element_text(
      angle = 52,
      hjust = 1,
      vjust = 1,
      face = "italic"
    ),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank(),
    axis.line.x.top = element_line(color = "black", linewidth = 0.4),
    axis.line.y.left = element_line(color = "black", linewidth = 0.4),
    axis.ticks.y.left = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"
  ) +
  facet_wrap(~Species, nrow = 1, scales = "free_x")

# Custom legend row
year_leg <- tibble(
  x = c(1.3, 2.4),
  y = 1,
  fill = c("1989", "2023"),
  lab = c("1989", "2023")
)

gf_leg <- tibble(
  x = c(9.3, 11.3, 13),
  y = 1,
  image = c(gram, cush, forb),
  lab = c("Graminoid", "Cushion", "Other forb")
)

p_leg <- ggplot() +
  geom_tile(
    data = year_leg,
    aes(x = x, y = y, fill = fill),
    width = 0.3,
    height = 0.20,
    show.legend = FALSE
  ) +
  geom_text(
    data = year_leg,
    aes(x = x + 0.2, y = y, label = lab),
    hjust = 0,
    size = 4.5,
    family = "serif"
  ) +
  geom_image(
    data = gf_leg,
    aes(x = x, y = y, image = image),
    inherit.aes = FALSE,
    size = 0.5,
    by = "width"
  ) +
  geom_text(
    data = gf_leg,
    aes(x = x + 0.3, y = y, label = lab),
    hjust = 0,
    size = 4.5,
    family = "serif"
  ) +
  scale_fill_manual(values = c("1989" = "skyblue3", "2023" = "red4")) +
  coord_cartesian(xlim = c(0.2, 14.3), ylim = c(0.82, 1.18), clip = "off") +
  theme_void()

fig1 <- p_leg / p_main + plot_layout(heights = c(0.12, 1))

ggsave("figures/barplot.png", fig1, width = 10, height = 5, dpi = 600)
