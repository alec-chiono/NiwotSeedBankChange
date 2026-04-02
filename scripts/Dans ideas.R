
librarian::shelf(tidyverse, ggbreak, ggthemes, ggimage, patchwork, ggdist)

# DATA -------------------------------------------------------------------------
## Download
source("scripts/source/download_data.R")

plot_df <- data_list$seedbank_composition.ac_hh.data %>%
  filter(
    substr(USDA_code, 1, 1) != 2 & USDA_code != "CAREX" & USDA_code != "POA"
  ) %>%
  mutate(Species = sub("\\s*var\\..*$", "", USDA_name), year = factor(year)) %>%
  group_by(year, Species) %>%
  summarize(sum = sum(count), .groups = "drop")

species_order <- plot_df %>%
  filter(year == "1989") %>%
  arrange(desc(sum)) %>%
  pull(Species)

# 1. Build lookup for symbol file paths
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

icon_df <- tibble(
  Species = names(symbol_map),
  image = unname(symbol_map),
  y_icon = -2.2
) %>%
  mutate(Species = factor(Species, levels = species_order))

f1 <- plot_df %>%
  mutate(Species = factor(Species, levels = species_order)) %>%
  ggplot(aes(x = Species, y = sum, color = year, fill = year)) +
  geom_col(position = "dodge") +
  geom_image(
    data = icon_df,
    aes(x = Species, y = y_icon, image = image),
    inherit.aes = FALSE,
    size = 0.065,
    by = "width"
  ) +
  scale_y_continuous(
    name = "# Germinable Seeds",
    expand = expansion(mult = c(0, 0.05), add = c(3.2, 0))
  ) +
  scale_fill_manual(
    name = "Year",
    values = c("skyblue3", "red4"),
    aesthetics = c("fill", "color")
  ) +
  scale_y_break(c(21, 84), scales = 0.25) +
  coord_cartesian(clip = "off") +
  theme_tufte() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      face = "italic"
    ),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.4),
    axis.line.y.left = element_line(color = "black", linewidth = 0.4),
    axis.ticks.y.left = element_line(color = "black"),
    strip.text = element_blank(),
    legend.position = "top",
    plot.margin = margin(5.5, 5.5, 55, 5.5)
  ) +
  facet_wrap(~Species, nrow = 1, scales = "free_x")
ggsave("figures/barplot.png", f1, width = 10, height = 5, dpi = 600)

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
    expand = expansion(mult = c(0, 0.05), add = c(3.2, 0))
  ) +
  #scale_x_discrete(expand = expansion(add = c(2.5, 0.2))) +
  scale_y_break(c(21, 84), scales = 0.25) +
  coord_cartesian(clip = "off") +
  theme_tufte() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 0.84,
      vjust = 0.93,
      face = "italic"
    ),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.4),
    axis.line.y.left = element_line(color = "black", linewidth = 0.4),
    axis.ticks.y.left = element_line(color = "black"),
    strip.text = element_blank(),
    legend.position = "none",
    plot.margin = margin(5.5, 5.5, 55, 5.5)
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
  annotate(
    "text",
    x = 0.35,
    y = 1,
    label = "Year",
    hjust = 0,
    size = 6,
    family = "serif"
  ) +
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
  annotate(
    "text",
    x = 6.7,
    y = 1,
    label = "Growth form",
    hjust = 0,
    size = 6,
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

f1 <- p_leg / p_main + plot_layout(heights = c(0.12, 1))

ggsave("figures/barplot.png", f1, width = 10, height = 5, dpi = 600)

## Plot by habitat
hab_df <- data_list$seedbank_composition.ac_hh.data %>%
  filter(
    substr(USDA_code, 1, 1) != 2 & USDA_code != "CAREX" & USDA_code != "POA"
  ) %>%
  mutate(
    Species = sub("\\s*var\\..*$", "", USDA_name),
    year = factor(year)
  ) %>%
  group_by(year, habitat, Species) %>%
  summarize(sum = sum(count), .groups = "drop")

make_hab_plot <- function(dat, hab, show_x = TRUE) {
  hab_order <- dat %>%
    filter(habitat == hab, year == "1989") %>%
    arrange(desc(sum)) %>%
    pull(Species)

  dat %>%
    filter(habitat == hab) %>%
    mutate(Species = factor(Species, levels = hab_order)) %>%
    ggplot(aes(x = Species, y = sum, fill = year, color = year)) +
    geom_col(position = "dodge") +
    coord_cartesian(ylim = c(0, 25)) +
    scale_y_continuous(name = "# Germinable Seeds") +
    scale_fill_manual(
      name = "Year",
      values = c("skyblue3", "red4"),
      aesthetics = c("fill", "color")
    ) +
    labs(title = hab, x = if (show_x) "Species" else NULL) +
    ggthemes::theme_tufte() +
    theme(
      plot.title = element_text(hjust = 0, face = "bold"),
      axis.text.x = if (show_x) {
        element_text(angle = 45, hjust = 1, face = "italic")
      } else {
        element_blank()
      },
      axis.ticks.x.bottom = if (show_x) {
        element_line(color = "black")
      } else {
        element_blank()
      },
      axis.line.x.bottom = element_line(color = "black", linewidth = 0.4),
      axis.line.y.left = element_line(color = "black", linewidth = 0.4),
      axis.ticks.y.left = element_line(color = "black"),
      legend.position = "top"
    )
}

p_mesic <- make_hab_plot(hab_df, "mesic", show_x = FALSE)
p_xeric <- make_hab_plot(hab_df, "xeric", show_x = TRUE)

(p_mesic / p_xeric) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

## Dans Fig 4 version
left_join(seed_draws, veg_draws, by=join_by(.draw, habitat, USDA_name)) %>%
  rename(seed_change=value.x, veg_change=value.y) %>%
  ggplot(aes(x=veg_change, y=seed_change, color=habitat)) +
  geom_hline(yintercept=0, linetype=2, color="grey") +
  geom_vline(xintercept=0, linetype=2, color="grey") +
  geom_point(alpha=.02) +
  scale_color_manual(values=c("mesic"="green4", "xeric"="tan4"), name="Habitat") +
  facet_wrap(~USDA_name, nrow=3, scales="free") +
  theme(strip.text = element_text(face = "italic"))

rel_draws %>%
  ggplot(aes(x=b_seed, fill=habitat)) +
  stat_slab(alpha=0.5, normalize="groups") +
  geom_vline(xintercept=0, linetype=2, color="red") +
  scale_x_continuous(name="Relationship between Seed Change and Veg Change") +
  scale_fill_manual(values=c("mesic"="green4", "xeric"="tan4"), guide="none") +
  theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
  ) +
  facet_wrap(~USDA_name, nrow=3) +
  theme(strip.text = element_text(face = "italic"))
