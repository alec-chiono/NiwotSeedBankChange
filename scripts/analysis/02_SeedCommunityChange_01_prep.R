# Prep data for analysis on seed bank community change
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ----
library(tidyverse)

# DATA ----
# Source data download function
source("scripts/source/download_data.R")

# Download and wrangle data
sc_df <-
  download_data(seed = TRUE)$seedbank_composition.ac_hh.data %>%
  filter(
    substr(USDA_code, 1, 1) != 2 & USDA_code != "CAREX" & USDA_code != "POA"
  ) %>% #remove records not identified to species
  select(year:depth, USDA_name, count) %>% #select relevant columns
  pivot_wider(names_from = USDA_name, values_from = count) %>% #pivot into proper format for analyses
  filter(rowSums(across(-(year:depth))) > 0) #remove samples that had no species present

# Pull out abundance data into separate matrix
sc_matrix <-
  sc_df %>%
  select(-(year:depth)) %>%
  as.matrix()

# Pull out meta data columns into separate data.frame
sc_meta <-
  sc_df %>%
  select(year:depth) %>%
  mutate(across(everything(), as.factor)) #pull out just metadata and make year a factor
