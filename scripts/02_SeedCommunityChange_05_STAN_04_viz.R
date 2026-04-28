# Visualizations of posteriors for Bayesian model for estimating richness, evenness, and diversity
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ----
library(tidybayes)
library(ggdist)
library(patchwork)

# MODEL ----
source("scripts/02_SeedCommunityChange_05_STAN_02_fit.R")

# POSTERIOR PREDICTIONS ----
# Retrieve posteriors for seed counts
post2_counts <-
  tidy_draws(fit2) %>%
  select(.draw, starts_with("predicted_seed_counts")) %>%
  pivot_longer(
    -.draw,
    names_to = c("metric", "year", "habitat", "species"),
    names_pattern = "(predicted_seed_counts)\\[(\\d+),(\\d+),(\\d+)\\]",
    values_to = "raw_count"
  ) %>%
  group_by(species) %>%
  mutate(scaled_count = raw_count / sd(raw_count)) %>%
  group_by(.draw, year, habitat) %>%
  summarize(scaled_total = sum(scaled_count), .groups = "drop") %>%
  mutate(
    metric = "Total Seeds",
    year = factor(
      year,
      levels = 1:dlist$N_years,
      labels = sort(unique(seedbank_df$year))
    ),
    habitat = factor(
      habitat,
      levels = 1:dlist$N_habitats,
      labels = sort(unique(seedbank_df$habitat))
    ),
    value = scaled_total
  ) %>%
  select(.draw, metric, year, habitat, value)

# Retrieve posteriors for richness, evenness, and diversity indices
post2_indices <-
  tidy_draws(fit2) %>%
  select(
    .draw,
    starts_with("richness"),
    starts_with("evenness"),
    starts_with("hill_N1")
  ) %>%
  pivot_longer(
    -.draw,
    names_to = c("metric", "year", "habitat"),
    names_pattern = "(predicted_seed_counts|richness|evenness|hill_N1)\\[(\\d+),(\\d+)]",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("richness", "evenness", "hill_N1"),
      labels = c(
        "Richness",
        "Evenness",
        "Diversity"
      )
    ),
    year = factor(
      year,
      levels = 1:dlist$N_years,
      labels = sort(unique(seedbank_df$year))
    ),
    habitat = factor(
      habitat,
      levels = 1:dlist$N_habitats,
      labels = sort(unique(seedbank_df$habitat))
    )
  )

# Collate posteriors together
post2 <- bind_rows(post2_indices, post2_counts) %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "Total Seeds",
        "Richness",
        "Evenness",
        "Diversity"
      )
    )
  )

# VIZ ----
# Set plotting theme
theme_set(
  ggthemes::theme_tufte() + theme(panel.border = element_rect(fill = NA))
)

## Fig. 1B: Estimates for each year and habitat
fig1B <- post2 %>%
  ggplot(aes(x = year)) +
  stat_slab(
    aes(y = value),
    fill = "black",
    color = "black",
    linewidth = 0.1,
    normalize = "groups"
  ) +
  facet_grid(metric ~ habitat, scales = "free_y") +
  scale_y_continuous(name = "Estimated Value", limits = c(0, NA))

## Fig. 1C Estimates for change in indices over time
fig1C <- post2 %>%
  arrange(.draw, metric, habitat, year) %>%
  group_by(.draw, metric, habitat) %>%
  mutate(diff = lead(value) - value, .groups = "drop") %>%
  filter(!is.na(diff)) %>%
  ggplot() +
  stat_slab(
    aes(y = diff),
    fill = "black",
    color = "black",
    linewidth = 0.1,
    normalize = "groups"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid(metric ~ habitat, scales = "free_y") +
  scale_y_continuous(name = "Difference between Years") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

## Get ggplot obj for Fig. 1A
suppressMessages(capture.output(
  source("scripts/analysis/02_SeedCommunityChange_04_RDA.R"),
  file = nullfile()
))

## Collate Figure 1
fig1 <- fig1A / (fig1B + fig1C) + plot_annotation(tag_levels = 'A')

## Write Figure 1
ggsave(
  "figures/fig1.pdf",
  fig1,
  width = 7.5,
  height = 10,
  units = "in",
  dpi = 600
)
