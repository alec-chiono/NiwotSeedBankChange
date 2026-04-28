# Compare change in seed bank community composition over time using PERMANOVA and BETADISPER
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ----
library(vegan)

# DATA ----
source("scripts/02_SeedCommunityChange_01_prep.R")

# Calculate Bray-Curtis distance matrix
bc_dist <- vegdist(sc_matrix, method = "bray")

# PERMANOVA ----
adonis2(
  bc_dist ~ year * habitat * depth,
  data = sc_meta,
  permutations = how(nperm = 999, blocks = sc_meta$plot),
  by = "terms"
)

# BETADISPER ----
## Test for homogeneity of multivariate dispersions
# By year
bd_year <- betadisper(bc_dist, group = sc_meta$year)
permutest(bd_year, permutations = 999)

# By habitat
bd_habitat <- betadisper(bc_dist, group = sc_meta$habitat)
permutest(bd_habitat, permutations = 999)

# By depth
bd_depth <- betadisper(bc_dist, group = sc_meta$depth)
permutest(bd_depth, permutations = 999)
