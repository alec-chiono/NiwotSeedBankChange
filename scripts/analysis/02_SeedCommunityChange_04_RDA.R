# Compare change in seed bank community composition over time using redundancy analysis (RDA)
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ----
library(vegan)

# DATA ----
source("scripts/analysis/02_SeedCommunityChange_01_prep.R")

# DCA ----
# to determine if we should do RDA or CCA
dca <- decorana(sc_matrix)
print(dca) # DCA1 axis length is 1.003, so RDA is appropriate

# RDA ----
sc_rda <- rda(
  decostand(sc_matrix, method = "hellinger") ~ #Hellinger transformation (https://r.qcbs.ca/workshop09/book-en/transformations.html)
    year * habitat * depth + Condition(plot),
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
) #check which terms are significant

# VIZ ----
# Set plotting theme
theme_set(
  ggthemes::theme_tufte() + theme(panel.border = element_rect(fill = NA))
)

# Extract sample scores
sc_scores <- sc_meta %>%
  cbind(scores(sc_rda)$sites)

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
