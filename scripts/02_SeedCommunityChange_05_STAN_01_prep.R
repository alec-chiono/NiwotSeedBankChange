# Prep data for Bayesian model for estimating richness, evenness, and diversity
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ----
library(tidyverse)

# DATA ----
# Source data download function
source("scripts/00_DownloadData.R")

# Download and wrangle data
seedbank_df <- download_data(seed = TRUE)$seedbank_composition.ac_hh.data %>%
  select(year:USDA_code, count) %>% #select relevant columns
  filter(
    substr(USDA_code, 1, 1) != 2 & USDA_code != "CAREX" & USDA_code != "POA"
  ) %>% #remove records not identified to species
  group_by(year, habitat, plot, USDA_code) %>%
  summarize(count = sum(count), .groups = "drop")

# Put data into list for cmdstanr model
dlist <- list(
  N_obs = nrow(seedbank_df),
  N_species = length(unique(seedbank_df$USDA_code)),
  N_plots = length(unique(seedbank_df$plot)),
  N_years = length(unique(seedbank_df$year)),
  N_habitats = length(unique(seedbank_df$habitat)),
  count = seedbank_df$count,
  species = as.integer(factor(seedbank_df$USDA_code)),
  plot = as.integer(factor(seedbank_df$plot)),
  year_idx = as.integer(factor(seedbank_df$year)),
  habitat_idx = as.integer(factor(seedbank_df$habitat))
)
