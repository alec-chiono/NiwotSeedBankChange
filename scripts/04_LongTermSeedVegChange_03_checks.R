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

# POSTERIOR RETRODICTIVE CHECKS ----
source("scripts/00_Function_plot_prc_hist.R")

# Figure S2A: are the veg slopes reasonable?
figS2A <- plot_bveg_vs_ols(fit4, vegf, USDA_lookup)

# Figure S2B: does the seed predictive distribution differ by veg trend?
figS2B <- plot_seed_prc(fit4, seedf)

# Collate Figure S2
figS2 <- figS2A / figS2B + plot_annotation(tag_levels = "A")

# Write out Figure S2
ggsave(
  "figures/FigureS2.pdf",
  figS2,
  width = 6.5,
  height = 6.5,
  units = "in",
  dpi = 600
)
