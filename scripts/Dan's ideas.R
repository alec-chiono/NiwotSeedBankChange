
librarian::shelf(tidyverse, ggbreak, ggthemes, patchwork)

# DATA -------------------------------------------------------------------------
## Download
source("scripts/source/download_data.R")

plot_df <- data_list$seedbank_composition.ac_hh.data %>%
  filter(substr(USDA_code, 1, 1) != 2 & USDA_code != "CAREX" & USDA_code != "POA") %>%
  mutate(Species = sub("\\s*var\\..*$", "", USDA_name),
         year = factor(year)) %>%
  group_by(year, Species) %>%
  summarize(sum = sum(count), .groups = "drop")

species_order <- plot_df %>%
  filter(year == "1989") %>%
  arrange(desc(sum)) %>%
  pull(Species)

plot_df %>%
  mutate(Species = factor(Species, levels = species_order)) %>%
  ggplot(aes(x = Species, y = sum, color = year, fill = year)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(name = "# Germinable Seeds") +
  scale_fill_manual(name = "Year", values = c("skyblue3", "red4"),
                    aesthetics = c("fill", "color")) +
  scale_y_break(c(21, 84), scales = 0.25) +
  ggthemes::theme_tufte() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.4),
    axis.line.y.left   = element_line(color = "black", linewidth = 0.4),
    axis.ticks.y.left  = element_line(color = "black"),
    strip.text = element_blank(),
    legend.position = "top"
  ) +
  facet_wrap(~Species, nrow = 1, scales = "free_x")

## Plot by habitat
hab_df <- data_list$seedbank_composition.ac_hh.data %>%
  filter(substr(USDA_code, 1, 1) != 2 & USDA_code != "CAREX" & USDA_code != "POA") %>%
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
      axis.text.x = if (show_x) element_text(angle = 45, hjust = 1, face = "italic") else element_blank(),
      axis.ticks.x.bottom = if (show_x) element_line(color = "black") else element_blank(),
      axis.line.x.bottom = element_line(color = "black", linewidth = 0.4),
      axis.line.y.left   = element_line(color = "black", linewidth = 0.4),
      axis.ticks.y.left  = element_line(color = "black"),
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
  mutate(Species=sub("\\s*var\\..*$", "", USDA_name)) %>%
  ggplot(aes(x=veg_change, y=seed_change, color=habitat)) +
  geom_hline(yintercept=0, linetype=2, color="grey") +
  geom_vline(xintercept=0, linetype=2, color="grey") +
  geom_point(alpha=.1) +
  scale_color_manual(values=c("mesic"="green4", "xeric"="tan4"), name="Habitat") +
  facet_wrap(~Species, nrow=3) +
  theme(strip.text = element_text(face = "italic"))
