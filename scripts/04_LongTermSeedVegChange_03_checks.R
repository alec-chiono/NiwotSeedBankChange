# Diagnostics and posterior predictive checks for Bayesian model evaluating relationship between long-term above-ground vegetation changes and seed bank changes
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ----
library(bayesplot)
library(posterior)
library(patchwork)

# MODEL ----
source("scripts/04_LongTermSeedVegChange_02_fit.R")

# DIAGNOSTICS ----
## cmdstanr doesn't automatically warn you when there are Rhat issues
fit4$cmdstan_diagnose()

# RANK PLOTS ----
# Find parameters that had poor mixing
summ <- fit4$draws() %>% summarise_draws()
poor_mixing <- summ %>%
  filter(rhat > 1.01 | ess_bulk < 400 | ess_tail < 400) %>%
  arrange(desc(rhat)) %>%
  pull(variable)

# Rank plots for poorly mixed parameters, if any
if (length(poor_mixing) == 0) {
  message("No parameters with poor mixing detected.")
} else {
  mcmc_rank_overlay(
    fit4$draws(),
    pars = poor_mixing
  )
}


# PPC ----
## Posterior Predictive Check
### Extract draws
seed_rep <- fit4$draws("seed_change_rep", format = "matrix")
veg_rep <- fit4$draws("veg_count_rep", format = "matrix")
y_seed <- dlist$seed_change
y_veg <- dlist$veg_count

### Grouped by habitat and species
ppc4A <- ppc_stat_grouped(
  y_seed,
  seed_rep,
  group = USDA_lookup$USDA_code[dlist$seed_species],
  stat = "mean"
)
ppc4B <- ppc_stat_grouped(
  y_seed,
  seed_rep,
  group = USDA_lookup$USDA_code[dlist$seed_species],
  stat = "sd"
)
ppc4C <- ppc_stat_grouped(
  y_veg,
  veg_rep,
  group = USDA_lookup$USDA_code[dlist$veg_species],
  stat = "mean"
)
ppc4D <- ppc_stat_grouped(
  y_veg,
  veg_rep,
  group = USDA_lookup$USDA_code[dlist$veg_species],
  stat = "sd"
)

ppc4 <- (ppc4A + ppc4B) / (ppc4C + ppc4D) + plot_annotation(tag_level = "A")
ggsave(
  "figures/figS2.pdf",
  ppc4,
  width = 10,
  height = 10,
  units = "in",
  dpi = 600
)
