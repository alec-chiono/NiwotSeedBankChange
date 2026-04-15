library(tidyverse, ggbreak, ggthemes, ggimage, patchwork, ggdist)

# DATA -------------------------------------------------------------------------
## Download
source("scripts/source/download_data.R")

gram <- "images/graminoid.png"
cush <- "images/cushion.png"
forb <- "images/forb.png"

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

# icons under species names
icon_df <- tibble(
  Species = names(symbol_map),
  image = unname(symbol_map),
  y_icon = -2.2
) %>%
  mutate(Species = factor(Species, levels = species_order))

# main plot
p_main <- plot_df %>%
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

# custom legend row
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
  # annotate(
  #   "text",
  #   x = 0.35,
  #   y = 1,
  #   label = "Year",
  #   hjust = 0,
  #   size = 6,
  #   family = "serif"
  # ) +
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
  # annotate(
  #   "text",
  #   x = 6.7,
  #   y = 1,
  #   label = "Growth form",
  #   hjust = 0,
  #   size = 6,
  #   family = "serif"
  # ) +
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

f1 <- p_leg / p_main + plot_layout(heights = c(0.12, 1))

ggsave("figures/barplot.png", f1, width = 10, height = 5, dpi = 600)
