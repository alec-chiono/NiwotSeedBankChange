# Prep data for analysis comparing co-located above-ground vegetation and seed bank community composition
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ----
library(tidyverse)

# DATA ----
# Source data download function
source("scripts/00_Function_data_download.R")

# Dowload data
data_list <- download_data(seed = TRUE, coloc_veg = TRUE)

### Seed bank data
seed_df <- data_list$seedbank_composition.ac_hh.data %>%
  filter(
    year == 2023, #use only recent data
    substr(USDA_code, 1, 1) != 2 & USDA_code != "CAREX" & USDA_code != "POA" #remove records not identified to species
  ) %>%
  group_by(habitat, plot, replicate, USDA_name) %>%
  summarize(count = sum(count), .groups = "drop") %>% #sum counts across depths
  mutate(
    community = "seed", #denote this as seed community data
    present = ifelse(count > 0, 1, 0) #convert count data into presence-absence data
  ) %>%
  select(community, habitat:USDA_name, present) #select relevant columns

### Co-located veg data
veg_df <- data_list$veg_composition.ac_hh.data %>%
  filter(
    substr(USDA_code, 1, 1) != 2 & USDA_code != "POA" & USDA_code != "CAREX"
  ) %>% #remove records not identified to species
  mutate(community = "vegetation") %>% #denote this as vegetation community data
  select(community, habitat:plot, USDA_name) %>% #select relevant columns
  distinct() %>%
  mutate(present = 1) #convert into presence data

### Combine seed and veg data
sv_df <- bind_rows(seed_df, veg_df) %>%
  complete(
    nesting(community, habitat, plot, replicate),
    USDA_name,
    fill = list(present = 0)
  ) %>% #turn into complete presence-absence data set
  pivot_wider(names_from = USDA_name, values_from = present) %>% #pivot into wide format for use in RDA
  mutate(across(community:replicate, as.factor)) #make metadata columns factors
