# Create Table S1, summary table of raw seed and veg data
# Alec Chiono; alec.chiono@colorado.edu

librarian::shelf(tidyverse, flextable, officer) #load packages

# Load data
source("scripts/source/download_data.R")

# ── Table S1: Seedbank ──────────────────────────────────────────────────────

seedbank_wide <- data_list$seedbank_composition.ac_hh.data %>%
  filter(
    substr(USDA_code, 1, 1) != 2 & USDA_code != "CAREX" & USDA_code != "POA"
  ) %>%
  group_by(USDA_name, habitat, year, depth) %>%
  summarize(
    n = n(),
    mean_count = mean(count),
    sd_count = sd(count),
    sem = sd_count / sqrt(n),
    mean_sem = paste0(round(mean_count, 1), " \u00b1 ", round(sem, 1)),
    .groups = "drop"
  ) %>%
  unite("group_key", c(habitat, year, depth), sep = "_") %>%
  select(USDA_name, group_key, mean_sem) %>%
  pivot_wider(
    names_from = group_key,
    values_from = mean_sem,
    values_fill = "\u2013"
  )

# Rename columns to display labels
col_rename <- setNames(
  names(seedbank_wide),
  names(seedbank_wide) %>%
    str_replace(".*0to5.*", "0–5 cm") %>%
    str_replace(".*5to10.*", "5–10 cm")
)

ft_s1 <- flextable(seedbank_wide) %>%
  set_header_labels(
    values = setNames(
      # values = new display labels
      names(seedbank_wide) %>%
        str_replace("^(mesic|xeric)_(1989|2023)_0to5$", "0\u20135 cm") %>%
        str_replace("^(mesic|xeric)_(1989|2023)_5to10$", "5\u201310 cm") %>%
        str_replace("^USDA_name$", "Species"),
      # names = original column names (used as keys)
      names(seedbank_wide)
    )
  ) %>%
  italic(j = "USDA_name") %>%
  # Level 2 spanners: Mesic / Xeric
  add_header_row(
    values = c("", "1989", "2023", "1989", "2023"),
    colwidths = c(
      1,
      sum(
        grepl("mesic", names(seedbank_wide)[-1]) &
          grepl("1989", names(seedbank_wide)[-1])
      ),
      sum(
        grepl("mesic", names(seedbank_wide)[-1]) &
          grepl("2023", names(seedbank_wide)[-1])
      ),
      sum(
        grepl("xeric", names(seedbank_wide)[-1]) &
          grepl("1989", names(seedbank_wide)[-1])
      ),
      sum(
        grepl("xeric", names(seedbank_wide)[-1]) &
          grepl("2023", names(seedbank_wide)[-1])
      )
    )
  ) %>%
  # Level 3 spanners: 1989 / 2023 within each habitat
  add_header_row(
    values = c("", "Mesic", "", "Xeric", ""),
    colwidths = c(
      1,
      sum(grepl("mesic_1989", names(seedbank_wide)[-1])),
      sum(grepl("mesic_2023", names(seedbank_wide)[-1])),
      sum(grepl("xeric_1989", names(seedbank_wide)[-1])),
      sum(grepl("xeric_2023", names(seedbank_wide)[-1]))
    )
  ) %>%
  align(align = "center", part = "all") %>%
  align(j = "USDA_name", align = "left", part = "all") %>%
  bold(part = "header") %>%
  set_caption(
    "Table S1. Mean \u00b1 SE of seed density per four 221 cm\u00b3 cores"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  fit_to_width(max_width = 8)

save_as_docx(ft_s1, path = "tables/tableS1.docx")


# ── Table S2: Vegetation ────────────────────────────────────────────────────

veg_wide <- data_list$veg_composition.ac_hh.data %>%
  filter(
    substr(USDA_code, 1, 1) != 2 & USDA_code != "CAREX" & USDA_code != "POA"
  ) %>%
  mutate(hit_value = if_else(hit_type == "extra", 0.5, 1)) %>%
  group_by(USDA_name, habitat, plot) %>%
  summarize(count = sum(hit_value), .groups = "drop") %>%
  complete(USDA_name, plot, habitat, fill = list(count = 0)) %>%
  group_by(USDA_name, habitat) %>%
  summarize(
    n = n(),
    mean_count = mean(count),
    sd_count = sd(count),
    sem = sd_count / sqrt(n),
    mean_sem = paste0(round(mean_count, 1), " \u00b1 ", round(sem, 1)),
    .groups = "drop"
  ) %>%
  select(USDA_name, habitat, mean_sem) %>%
  pivot_wider(names_from = habitat, values_from = mean_sem)

ft_s2 <- flextable(veg_wide) %>%
  set_header_labels(
    USDA_name = "Species",
    mesic = "Mesic",
    xeric = "Xeric"
  ) %>%
  italic(j = "USDA_name") %>%
  align(align = "center", part = "all") %>%
  align(j = "USDA_name", align = "left", part = "all") %>%
  bold(part = "header") %>%
  set_caption(
    "Table S2. Mean \u00b1 SE of vegetation hits per 1 m\u00b2 plot in 2024."
  ) %>%
  theme_booktabs() %>%
  autofit()

save_as_docx(ft_s2, path = "tables/tableS2.docx")
