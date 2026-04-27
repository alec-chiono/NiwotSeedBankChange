# Compare change in seed bank community composition over time using redundancy analysis (RDA)
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ---------------------------------------------------------------------
library(tidyverse)
library(vegan)
library(betapart)

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

# --- Helper: extract mean dissimilarity between two groups ---
mean_between <- function(dist_obj, meta, var, group1, group2) {
  mat <- as.matrix(dist_obj)
  idx1 <- which(meta[[var]] == group1)
  idx2 <- which(meta[[var]] == group2)
  mean(mat[idx1, idx2])
}

# --- Define your four comparisons ---
# Adjust variable names and levels to match your actual sc_meta columns

# Subset metadata indices for each year
meta_1989 <- sc_meta %>% mutate(.idx = row_number()) %>% filter(year == 1989)
meta_2023 <- sc_meta %>% mutate(.idx = row_number()) %>% filter(year == 2023)

# Helper for subsetting a dist to specific indices then taking between-group mean
mean_between_sub <- function(dist_obj, meta_sub, var, g1, g2) {
  mat <- as.matrix(dist_obj)[meta_sub$.idx, meta_sub$.idx]
  idx1 <- which(meta_sub[[var]] == g1)
  idx2 <- which(meta_sub[[var]] == g2)
  mean(mat[idx1, idx2])
}

comparisons <- tribble(
  ~group    , ~label                    , ~year_filter ,
  "Habitat" , "Mesic vs. xeric\n(1989)" ,         1989 ,
  "Habitat" , "Mesic vs. xeric\n(2023)" ,         2023 ,
  "Year"    , "1989 vs. 2023\n(xeric)"  , NA           ,
  "Year"    , "1989 vs. 2023\n(mesic)"  , NA
)

# Build results manually for each comparison
results <- bind_rows(
  # Habitat comparisons (within year)
  tibble(
    group = "Habitat",
    label = "Mesic vs. xeric\n(1989)",
    component = c("Species replacement", "Abundance gradient"),
    value = c(
      mean_between_sub(
        beta_abund$beta.bray.bal,
        sc_meta %>% mutate(.idx = row_number()) %>% filter(year == 1989),
        "habitat",
        "mesic",
        "xeric"
      ),
      mean_between_sub(
        beta_abund$beta.bray.gra,
        sc_meta %>% mutate(.idx = row_number()) %>% filter(year == 1989),
        "habitat",
        "mesic",
        "xeric"
      )
    )
  ),
  tibble(
    group = "Habitat",
    label = "Mesic vs. xeric\n(2023)",
    component = c("Species replacement", "Abundance gradient"),
    value = c(
      mean_between_sub(
        beta_abund$beta.bray.bal,
        sc_meta %>% mutate(.idx = row_number()) %>% filter(year == 2023),
        "habitat",
        "mesic",
        "xeric"
      ),
      mean_between_sub(
        beta_abund$beta.bray.gra,
        sc_meta %>% mutate(.idx = row_number()) %>% filter(year == 2023),
        "habitat",
        "mesic",
        "xeric"
      )
    )
  ),
  # Year comparisons (within habitat)
  tibble(
    group = "Year",
    label = "1989 vs. 2023\n(xeric)",
    component = c("Species replacement", "Abundance gradient"),
    value = c(
      mean_between(
        beta_abund$beta.bray.bal,
        sc_meta %>% filter(habitat == "xeric"),
        "year",
        1989,
        2023
      ),
      mean_between(
        beta_abund$beta.bray.gra,
        sc_meta %>% filter(habitat == "xeric"),
        "year",
        1989,
        2023
      )
    )
  ),
  tibble(
    group = "Year",
    label = "1989 vs. 2023\n(mesic)",
    component = c("Species replacement", "Abundance gradient"),
    value = c(
      mean_between(
        beta_abund$beta.bray.bal,
        sc_meta %>% filter(habitat == "mesic"),
        "year",
        1989,
        2023
      ),
      mean_between(
        beta_abund$beta.bray.gra,
        sc_meta %>% filter(habitat == "mesic"),
        "year",
        1989,
        2023
      )
    )
  )
) %>%
  mutate(
    group = factor(group, levels = c("Habitat", "Year")),
    label = factor(
      label,
      levels = c(
        "Mesic vs. xeric\n(1989)",
        "Mesic vs. xeric\n(2023)",
        "1989 vs. 2023\n(xeric)",
        "1989 vs. 2023\n(mesic)"
      )
    ),
    component = factor(
      component,
      levels = c("Species replacement", "Abundance gradient")
    )
  )

# --- Plot ---
figS1 <- ggplot(results, aes(x = label, y = value, fill = component)) +
  geom_col(width = 0.6) +
  facet_grid(~group, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = c(
      "Species replacement" = "#4E9AC7", # adjust to match your palette
      "Abundance gradient" = "#F4A460"
    ),
    name = NULL
  ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Mean pairwise dissimilarity (Bray–Curtis)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 10)
  )

ggsave(
  "figures/figS1.pdf",
  figS1,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)
