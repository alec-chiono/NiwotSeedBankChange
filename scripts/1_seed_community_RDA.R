# Compare change in seed bank community composition over time using redundancy analysis (RDA)
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ---------------------------------------------------------------------
librarian::shelf(tidyverse, vegan, betapart, indicspecies, flextable, officer)

# DATA -------------------------------------------------------------------------
## Download
source("scripts/source/download_data.R")

## Wrangle
sc_df <- data_list$seedbank_composition.ac_hh.data %>%
  filter(
    substr(USDA_code, 1, 1) != 2 & USDA_code != "CAREX" & USDA_code != "POA"
  ) %>% #remove records not identified to species
  select(year:depth, USDA_name, count) %>% #select relevant columns
  pivot_wider(names_from = USDA_name, values_from = count) %>% #pivot into proper format for analyses
  filter(rowSums(across(-(year:depth))) > 0) #remove samples that had no species present

# RDA ---------------------------------------------------------------------
## Pull out abundance data
sc_matrix <- sc_df %>%
  select(-(year:depth)) %>%
  as.matrix()

## Pull out meta data
sc_meta <- sc_df %>%
  select(year:depth) %>%
  mutate(across(everything(), as.factor)) #pull out just metadata and make year a factor

## DCA to determine if we should do RDA or CCA
dca <- decorana(sc_matrix)
print(dca) #DCA1 axis length is 1.003, so will do RDA

## Fit RDA
sc_rda <- rda(
  decostand(sc_matrix, method = "hellinger") ~ #Hellinger transformation (https://r.qcbs.ca/workshop09/book-en/transformations.html)
    year *
    habitat *
    depth + #variables of interest
    Condition(plot), #spatial structure
  data = sc_meta
)

## Evaluate model significance
anova.cca(sc_rda, step = 1000) #check overall significance
anova.cca(sc_rda, step = 1000, by = "axis") #check axes significance
anova.cca(
  sc_rda,
  step = 1000,
  by = "margin",
  scope = c(
    "year",
    "habitat",
    "depth",
    "year:habitat",
    "year:depth",
    "habitat:depth",
    "year:habitat:depth"
  )
) #check which terms significant

## Extract sample scores
sc_scores <- sc_meta %>%
  cbind(scores(sc_rda)$sites)


# FIGURES ----------------------------------------------------------------------
## Set plotting theme
theme_set(
  ggthemes::theme_tufte() + theme(panel.border = element_rect(fill = NA))
)

## Fig. 1A: Just sample scores
fig1A <- sc_scores %>%
  mutate(habitat_depth = paste(habitat, depth, sep = "_")) %>%
  ggplot(aes(x = RDA1, y = RDA2)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_point(aes(color = year, shape = habitat), size = 5) +
  scale_shape_manual(values = c(16, 17)) +
  scale_color_manual(values = c("skyblue3", "red4")) +
  theme(legend.position = "top")

ggsave(
  "figures/fig1A.png",
  fig1A,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)

# PERMANOVA --------------------------------------------------------------------
## Bray-Curtis distance matrix
bc_dist <- vegdist(sc_matrix, method = "bray")

## Fit PERMANOVA
adonis2(
  bc_dist ~
    year * habitat * depth,
  data = sc_meta,
  permutations = how(nperm = 999, blocks = sc_meta$plot),
  by = "terms"
)

# BETADISPER -------------------------------------------------------------------
## Test for homogeneity of multivariate dispersions (i.e., beta diversity)

# Year
bd_year <- betadisper(bc_dist, group = sc_meta$year)
permutest(bd_year, permutations = 999)

# Habitat
bd_habitat <- betadisper(bc_dist, group = sc_meta$habitat)
permutest(bd_habitat, permutations = 999)

# Depth
bd_depth <- betadisper(bc_dist, group = sc_meta$depth)
permutest(bd_depth, permutations = 999)

# BETAPART DECOMPOSITION -------------------------------------------------------
## Requires raw abundance matrix (not hellinger transformed)
beta_abund <- beta.pair.abund(sc_matrix, index.family = "bray")

## Total
adonis2(
  beta_abund$beta.bray ~ year * habitat * depth,
  data = sc_meta,
  permutations = how(nperm = 999, blocks = sc_meta$plot),
  by = "terms"
)

## Turnover component
adonis2(
  beta_abund$beta.bray.bal ~ year * habitat * depth,
  data = sc_meta,
  permutations = how(nperm = 999, blocks = sc_meta$plot),
  by = "terms"
)

## Abundance gradient component
adonis2(
  beta_abund$beta.bray.gra ~ year * habitat * depth,
  data = sc_meta,
  permutations = how(nperm = 999, blocks = sc_meta$plot),
  by = "terms"
)

## Plot results
### Index sets
idx <- list(
  "1989 vs. 2023\n(xeric)" = list(
    which(sc_meta$year == "1989" & sc_meta$habitat == "xeric"),
    which(sc_meta$year == "2023" & sc_meta$habitat == "xeric")
  ),
  "1989 vs. 2023\n(mesic)" = list(
    which(sc_meta$year == "1989" & sc_meta$habitat == "mesic"),
    which(sc_meta$year == "2023" & sc_meta$habitat == "mesic")
  ),
  "Mesic vs. xeric\n(1989)" = list(
    which(sc_meta$year == "1989" & sc_meta$habitat == "mesic"),
    which(sc_meta$year == "1989" & sc_meta$habitat == "xeric")
  ),
  "Mesic vs. xeric\n(2023)" = list(
    which(sc_meta$year == "2023" & sc_meta$habitat == "mesic"),
    which(sc_meta$year == "2023" & sc_meta$habitat == "xeric")
  )
)

### Extract matrices
bal_mat <- as.matrix(beta_abund$beta.bray.bal)
gra_mat <- as.matrix(beta_abund$beta.bray.gra)

# Build plot data
plot_data <- imap_dfr(
  idx,
  ~ tibble(
    comparison = .y,
    `Species replacement` = mean(bal_mat[.x[[1]], .x[[2]]]),
    `Abundance gradient` = mean(gra_mat[.x[[1]], .x[[2]]])
  )
) %>%
  pivot_longer(
    -comparison,
    names_to = "component",
    values_to = "dissimilarity"
  ) %>%
  mutate(
    component = factor(
      component,
      levels = c("Species replacement", "Abundance gradient")
    ),
    comparison = factor(comparison, levels = names(idx)),
    comparison_type = ifelse(
      grepl("vs\\. 2023|vs\\. 1989", comparison),
      "Year",
      "Habitat"
    )
  )

# Plot
figS1 <- ggplot(
  plot_data,
  aes(x = comparison, y = dissimilarity, fill = component)
) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c(
      "Species replacement" = "gray30",
      "Abundance gradient" = "gray75"
    ),
    name = NULL
  ) +
  facet_grid(~comparison_type, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Mean pairwise dissimilarity (Bray-Curtis)") +
  theme(legend.position = "top")

ggsave("figures/figS1.pdf", figS1, width = 6, height = 4.5, dpi = 300)

# INDICATOR SPECIES ANALYSIS ---------------------------------------------------
## Create a combined grouping variable for year x habitat combinations
sc_meta <- sc_meta %>%
  mutate(year_habitat = interaction(year, habitat, sep = "_"))

## Run indicator species analysis
## Using the abundance matrix (not presence/absence)
indval_result <- multipatt(
  sc_matrix,
  cluster = sc_meta$year_habitat,
  func = "IndVal.g", # corrected IndVal
  control = how(nperm = 999, blocks = sc_meta$plot)
)

summary(indval_result, indvalcomp = TRUE, alpha = 0.05)

## Write results to formatted tables
### Format p-values
indval_formatted <- indval_result$sign %>%
  rownames_to_column("species_code") %>%
  filter(p.value <= 0.05) %>%
  mutate(
    # Identify which group(s) the species is associated with
    group = case_when(
      `s.1989_mesic` == 1 &
        `s.2023_mesic` == 0 &
        `s.1989_xeric` == 0 &
        `s.2023_xeric` == 0 ~ "1989 — mesic only",
      `s.1989_xeric` == 1 &
        `s.2023_xeric` == 0 &
        `s.1989_mesic` == 0 &
        `s.2023_mesic` == 0 ~ "1989 — xeric only",
      `s.1989_mesic` == 1 &
        `s.1989_xeric` == 1 &
        `s.2023_mesic` == 0 &
        `s.2023_xeric` == 0 ~ "1989 — both habitats",
      `s.2023_mesic` == 1 &
        `s.2023_xeric` == 0 &
        `s.1989_mesic` == 0 &
        `s.1989_xeric` == 0 ~ "2023 — mesic only",
      `s.2023_xeric` == 1 &
        `s.2023_mesic` == 0 &
        `s.1989_mesic` == 0 &
        `s.1989_xeric` == 0 ~ "2023 — xeric only",
      `s.2023_mesic` == 1 &
        `s.2023_xeric` == 1 &
        `s.1989_mesic` == 0 &
        `s.1989_xeric` == 0 ~ "2023 — both habitats",
      TRUE ~ "multiple groups"
    ),
    # Round numeric columns
    stat = round(stat, 3),
    p.value = round(p.value, 3)
  ) %>%
  # Join species names if you have a lookup table
  # left_join(species_lookup, by = "species_code") %>%
  select(species_code, group, stat, p.value) %>%
  arrange(group, desc(stat)) %>%
  rename(
    "Species code" = species_code,
    "Group" = group,
    "IndVal.g" = stat,
    "p" = p.value
  ) %>%
  mutate(p = ifelse(p <= 0.001, "< 0.001", as.character(p)))

### Build flextable
ft <- flextable(indval_formatted) %>%
  # Italic species names if you have binomial names in the table
  italic(j = "Species code") %>%
  # Bold header
  bold(part = "header") %>%
  # Add horizontal lines to separate groups
  hline(i = which(diff(as.numeric(factor(indval_formatted$Group))) != 0)) %>%
  # Column widths
  width(j = "Species code", width = 1.5) %>%
  width(j = "Group", width = 2.0) %>%
  width(j = "IndVal.g", width = 0.8) %>%
  width(j = "p", width = 0.6) %>%
  # Align numeric columns
  align(j = c("IndVal.g", "p"), align = "center", part = "all") %>%
  # Add caption
  set_caption(
    "Table S3. Indicator species analysis results. Species significantly associated
              (p ≤ 0.05) with each year–habitat combination based on IndVal.g (Cáceres & Legendre 2009).
              IndVal.g ranges from 0 to 1, where higher values indicate stronger association."
  ) %>%
  theme_booktabs()

### Save to Word
save_as_docx(ft, path = "tables/tableS3.docx")
