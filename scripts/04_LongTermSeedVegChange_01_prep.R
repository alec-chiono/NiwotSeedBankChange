# Prep data for Bayesian model evaluating relationship between long-term above-ground vegetation changes and seed bank changes
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ----
library(tidyverse)

# DATA ----
# Source data download function
source("scripts/00_Function_data_download.R")

# Download long-term seed and veg data
data_list <- download_data(seed = TRUE, longterm_veg = TRUE)

# Load categorizations of long-term veg plots into mesic or xerix dry meadow
saddlegrid_habitat <- read.csv("data/saddptqd_xericmesic_categorization.csv")

# Wrangle seed bank data
seedbank_df <- data_list$seedbank_composition.ac_hh.data %>%
  select(year:USDA_name, count) %>% #select relevant columns
  filter(
    substr(USDA_code, 1, 1) != 2 & USDA_code != "CAREX" & USDA_code != "POA"
  ) %>% #remove records not identified to species
  group_by(year, habitat, plot, USDA_code, USDA_name) %>%
  summarize(count = sum(count), .groups = "drop") %>%
  group_by(habitat, plot, USDA_code) %>%
  filter(any(count > 0)) %>%
  ungroup() %>%
  group_by(USDA_code, USDA_name) %>%
  mutate(count_scaled = scale(count)[, 1]) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_wider(
    names_from = year,
    values_from = count_scaled,
    names_prefix = "count"
  ) %>%
  mutate(change = count2023 - count1989)

# Wrangle long-term veg data
saddlegrid_df <- data_list$saddptqd.hh.data.csv %>%
  filter(
    year >= 1989 & year <= 2023, #same time frame as seed bank study
    plot %in% saddlegrid_habitat$plot,
    hit_type %in% c("top", "bottom") #middle and extra hits haven't been done across whole time period
  ) %>%
  group_by(year, plot, USDA_code) %>%
  summarize(count = length(hit_type), .groups = "drop") %>%
  complete(plot, USDA_code, year, fill = list(count = 0)) %>% #add records for species not seen in some years
  group_by(plot, USDA_code) %>%
  filter(any(count > 0)) %>%
  ungroup() %>%
  group_by(plot, USDA_code) %>%
  mutate(
    count_scaled = scale(count, center = FALSE)[, 1],
    count_scaled = ifelse(is.na(count_scaled), 0, count_scaled)
  ) %>%
  ungroup() %>%
  merge(saddlegrid_habitat, ., by = "plot") #add habitat information for later use

## Remove species-habitat combinations that weren't actually present in seed bank data
spp_to_include <- data_list$seedbank_composition.ac_hh.data %>%
  filter(
    substr(USDA_code, 1, 1) != 2 & USDA_code != "CAREX" & USDA_code != "POA",
    count > 0
  ) %>%
  distinct(USDA_code, habitat) %>%
  select(USDA_code, habitat) %>%
  mutate(name = "seed", value = 1) %>%
  rbind(
    saddlegrid_df %>%
      filter(
        count > 0,
        USDA_code %in% seedbank_df$USDA_code
      ) %>%
      distinct(USDA_code, habitat) %>%
      select(USDA_code, habitat) %>%
      mutate(name = "veg", value = 1)
  ) %>%
  pivot_wider() %>%
  drop_na()

### Filter seed bank data
seedf <- seedbank_df %>%
  filter(
    paste0(USDA_code, habitat) %in%
      with(spp_to_include, paste0(USDA_code, habitat)) #remove species-habitat combos that weren't present in data
  )

### Filter veg data to only include species-habiat combos now in seedf
vegf <- saddlegrid_df %>%
  filter(
    paste0(USDA_code, habitat) %in%
      with(spp_to_include, paste0(USDA_code, habitat))
  ) #remove species-habitat combos that weren't present in data)

## Create lookup data.frames for variables that will be scaled or turned into IDs for stan
year_lookup <- vegf %>%
  mutate(year_scaled = scale(year)[, 1]) %>%
  select(year, year_scaled) %>%
  distinct(year, year_scaled) %>%
  arrange(year)

USDA_lookup <- seedf %>%
  mutate(
    USDA_code_id = as.integer(factor(USDA_code)),
    Species = sub("\\s*var\\..*$", "", USDA_name),
  ) %>%
  select(USDA_code_id, USDA_code, USDA_name, Species) %>%
  distinct() %>%
  arrange(USDA_code_id)

veg_plot_lookup <- vegf %>%
  mutate(
    plot_id = as.integer(factor(plot)),
    habitat_id = habitat == "xeric"
  ) %>%
  select(plot_id, plot, habitat_id, habitat) %>%
  distinct() %>%
  arrange(plot_id)

seed_plot_lookup <- seedf %>%
  mutate(
    plot_id = as.integer(factor(plot)),
    habitat_id = habitat == "xeric"
  ) %>%
  select(plot_id, plot, habitat_id, habitat) %>%
  distinct() %>%
  arrange(plot_id)

# Put data into list for cmdstanr model
dlist <- list(
  # Vegetation data
  Nyear = nrow(vegf),
  veg_year = vegf$year - min(vegf$year),
  Kspecies = nrow(USDA_lookup),
  Nveg = nrow(vegf),
  veg_species = as.integer(factor(vegf$USDA_code)),
  veg_count = vegf$count_scaled,
  Nveg_plots = nrow(veg_plot_lookup),
  veg_plot = as.integer(factor(vegf$plot)),
  veg_habitat = veg_plot_lookup$habitat_id,

  # Seed data
  Nseed = nrow(seedf),
  seed_change = seedf$change / (2023 - 1989),
  seed_species = as.integer(factor(seedf$USDA_code)),
  Nseed_plots = nrow(seed_plot_lookup),
  seed_plot = as.integer(factor(seedf$plot)),
  seed_habitat = seed_plot_lookup$habitat_id
)

stopifnot(
  # No seed obs where both years were zero (change = 0 from structural absence)
  # Note: change = 0 can still occur legitimately (equal non-zero scaled counts)
  !any(seedf$count1989 == 0 & seedf$count2023 == 0),
  # Stan array lengths are consistent
  length(dlist$veg_species) == dlist$Nveg,
  length(dlist$seed_species) == dlist$Nseed,
  length(dlist$veg_habitat) == dlist$Nveg_plots,
  length(dlist$seed_habitat) == dlist$Nseed_plots,
  # Species IDs are contiguous (no gaps from filtering)
  max(dlist$veg_species) == dlist$Kspecies,
  max(dlist$seed_species) == dlist$Kspecies
)
