# Create Table S1, summary table of raw seed and veg data
# Alec Chiono; alec.chiono@colorado.edu

librarian::shelf(tidyverse, gt, gtExtras) #load packages

# Load data
source("Zenodo_archiving/scripts/source/download_data.R")
saddlegrid_habitat <- read.csv("data/saddlegrid_habitat_types.csv")

data_list$seedbank_composition.ac_hh.data %>%
  filter(substr(USDA_code,1, 1)!=2 & USDA_code!="CAREX" & USDA_code!="POA") %>%  #remove records not identified to species
  group_by(USDA_name, habitat, year, depth) %>%
  summarize(
    n = n(),
    mean_count = mean(count),
    sd_count = sd(count),
    sem = sd_count / sqrt(n),
    mean_sem = paste0(round(mean_count, 1), " ± ", round(sem, 1)),
    .groups = 'drop'
  ) %>%
  unite("group_key", c(habitat, year, depth), sep = "_") %>%
  select(USDA_name, group_key, mean_sem) %>%
  pivot_wider(
    names_from = group_key,
    values_from = mean_sem,
    names_sep = "_",
    values_fill = "–"
  ) %>%
  gt(rowname_col = "USDA_name") %>%
  tab_header(title = html("Mean ± SE of seed density per 221 cm<sup>3</sup> sample")) %>%
  # L2: Habitat (unique ids)
  tab_spanner(label = "Mesic", columns = contains("mesic"), id = "mesic", level = 2) %>%
  tab_spanner(label = "Xeric", columns = contains("xeric"), id = "xeric", level = 2) %>%
  # L1: Year (unique ids per habitat)
  tab_spanner(label = "1989", columns = contains("mesic_1989"), id = "mesic_1989", level = 1) %>%
  tab_spanner(label = "1989", columns = contains("xeric_1989"), id = "xeric_1989", level = 1) %>%
  tab_spanner(label = "2023", columns = contains("mesic_2023"), id = "mesic_2023", level = 1) %>%
  tab_spanner(label = "2023", columns = contains("xeric_2023"), id = "xeric_2023", level = 1) %>%
  # Italic row names (species)
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_stub()
  ) %>%
  # Labels for 1989/2023 (update as needed)
  cols_label(
    contains("0to5") ~ "0-5cm",
    contains("5to10") ~ "5-10cm"
  ) %>%
  cols_align(align = "center") %>%
  sub_missing(missing_text = "–") %>%
  gtsave("Zenodo_archiving/figures/Table_S1.html")

# Veg
data_list$veg_composition.ac_hh.data %>%
  filter(
    substr(USDA_code,1, 1)!=2 & USDA_code!="CAREX" & USDA_code!="POA", #remove records not identified to species
    ) %>%
  mutate(hit_value=if_else(hit_type=="extra", 0.5, 1)) %>%
  group_by(USDA_name, habitat, plot) %>%
  summarize(count=sum(hit_value), .groups="drop") %>%
  complete(USDA_name, plot, habitat, fill=list(count=0)) %>%
  group_by(USDA_name, habitat) %>%
  summarize(
    n = n(),
    mean_count = mean(count),
    sd_count = sd(count),
    sem = sd_count / sqrt(n),
    mean_sem = paste0(round(mean_count, 1), " ± ", round(sem, 1)),
    .groups = 'drop'
  ) %>%
  select(USDA_name, habitat, mean_sem) %>%
  pivot_wider(
    names_from = habitat,
    values_from = mean_sem
  ) %>%
  gt(rowname_col = "USDA_name") %>%
  tab_header(title = html("Mean ± SE of vegetation hits per 1m<sup>2</sup> plot in 2024")) %>%
  # Italic row names (species)
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_stub()
  ) %>%
  # Labels for 1989/2023 (update as needed)
  cols_label(
    contains("mesic") ~ "Mesic",
    contains("xeric") ~ "Xeric"
  ) %>%
  cols_align(align = "center") %>%
  sub_missing(missing_text = "–") %>%
  gtsave("Zenodo_archiving/figures/Table_S2.html")
