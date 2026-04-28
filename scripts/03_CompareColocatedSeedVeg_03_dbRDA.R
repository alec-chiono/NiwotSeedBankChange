# Compare co-located above-ground vegetation and seed bank community composition using distance-based RDA
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ----
library(vegan)

# DATA ----
source("scripts/03_CompareColocatedSeedVeg_01_prep.R")

# Pull out presence/absence data and convert into Jaccard distance matrix
sv_dist <- sv_df %>%
  select(-(community:replicate)) %>%
  as.matrix() %>%
  vegdist(method = "jaccard")

# Pull out meta data
sv_meta <- sv_df %>%
  select((community:replicate)) %>%
  mutate(across(everything(), as.factor))

# dbRDA ----
sv_rda <- dbrda(
  sv_dist ~ community * habitat + Condition(plot),
  data = sv_meta,
  distance = "jaccard"
)

# Evaluate model significance
anova.cca(sv_rda, step = 1000) #check overall significance
anova.cca(sv_rda, step = 1000, by = "axis") #check axes significance
anova.cca(
  sv_rda,
  step = 1000,
  by = "margin",
  scope = c("community", "habitat", "community:habitat")
) #check which terms are significant


# VIZ ----
# Set plotting theme
theme_set(
  ggthemes::theme_tufte() + theme(panel.border = element_rect(fill = NA))
)

# Extract scores
sv_scores <- sv_meta %>%
  cbind(scores(sv_rda)$sites)

# Fig. 3: dbRDA comparing seed bank and veg community from same locations
fig3 <- sv_scores %>%
  ggplot(aes(x = dbRDA1, y = dbRDA2)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_point(aes(color = community, shape = habitat), size = 5) +
  scale_shape_manual(values = c(16, 17)) +
  scale_color_manual(values = c("brown4", "green4")) +
  theme(legend.position = "top")

# Write Figure 2
ggsave(
  "figures/fig3.pdf",
  fig3,
  width = 7.5,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)
