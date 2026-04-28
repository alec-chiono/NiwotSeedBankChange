# Compare co-located above-ground vegetation and seed bank community composition using PERMANOVA and BETADISPER
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ----
library(vegan)

# DATA ----
source("scripts/03_CompareColocatedSeedVeg_01_prep.R")

# Aggregate presence/absence data for soil samples for each plot so blocking design does not invalidate PERMANOVA
sv_df_agg <- sv_df %>%
  group_by(community, habitat, plot) %>%
  summarize(across(where(is.numeric), max), .groups = "drop")

# Convert presence/absence matrix into Jaccard distance matrix
sv_dist_agg <- sv_df_agg %>%
  select(-(community:plot)) %>%
  as.matrix() %>%
  vegdist(method = "jaccard")

# Pull out meta data
sv_meta_agg <- sv_df_agg %>%
  select((community:plot)) %>%
  mutate(across(everything(), as.factor))

# PERMANOVA ----
adonis2(
  sv_dist_agg ~ community * habitat,
  data = sv_meta_agg,
  permutations = how(blocks = sv_meta_agg$plot, nperm = 999),
  by = "terms"
)

# BETADISPER ----
## Test for homogeneity of multivariate dispersions
# By community (seed or veg)
bd_community <- betadisper(sv_dist_agg, group = sv_meta_agg$community)
permutest(bd_community, permutations = 999)
# By habitat
bd_habitat <- betadisper(sv_dist_agg, group = sv_meta_agg$habitat)
permutest(bd_habitat, permutations = 999)
