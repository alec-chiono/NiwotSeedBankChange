# Evaluate contributions of Species turnover and Abdundance shift to seed bank community change
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ----
library(vegan)
library(betapart)

# DATA ----
source("scripts/02_SeedCommunityChange_01_prep.R")

# BETAPART DECOMPOSITION ----
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

## Abdundance shift component
adonis2(
  beta_abund$beta.bray.gra ~ year * habitat * depth,
  data = sc_meta,
  permutations = how(nperm = 999, blocks = sc_meta$plot),
  by = "terms"
)

# VIZ ----
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
  "Habitat" , "Mesic vs. Xeric\n(1989)" ,         1989 ,
  "Habitat" , "Mesic vs. Xeric\n(2023)" ,         2023 ,
  "Year"    , "1989 vs. 2023\n(xeric)"  , NA           ,
  "Year"    , "1989 vs. 2023\n(mesic)"  , NA
)

# Build results manually for each comparison
results <- bind_rows(
  # Habitat comparisons (within year)
  tibble(
    group = "Habitat (Mesic vs. Xeric)",
    label = "1989",
    component = c("Species turnover", "Abdundance shift"),
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
    group = "Habitat (Mesic vs. Xeric)",
    label = "2023",
    component = c("Species turnover", "Abdundance shift"),
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
    group = "Year (1989 vs. 2023)",
    label = "Xeric",
    component = c("Species turnover", "Abdundance shift"),
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
    group = "Year (1989 vs. 2023)",
    label = "Mesic",
    component = c("Species turnover", "Abdundance shift"),
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
    group = factor(
      group,
      levels = c("Habitat (Mesic vs. Xeric)", "Year (1989 vs. 2023)")
    ),
    label = factor(
      label,
      levels = c(
        "1989",
        "2023",
        "Mesic",
        "Xeric"
      )
    ),
    component = factor(
      component,
      levels = c("Species turnover", "Abdundance shift")
    )
  )

# --- Plot ---
figS4 <- ggplot(results, aes(x = label, y = value, fill = component)) +
  geom_col(width = 0.6) +
  facet_grid(~group, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = c(
      "Species turnover" = "#4E9AC7", # adjust to match your palette
      "Abdundance shift" = "#F4A460"
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
  "figures/FigureS4.pdf",
  figS4,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)
