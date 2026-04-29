# Diagnostics and posterior predictive checks for Bayesian model for estimating richness, evenness, and diversity
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ----
library(bayesplot)
library(posterior)
library(loo)
library(patchwork)

# MODEL ----
source("scripts/02_SeedCommunityChange_05_STAN_02_fit.R")

# DIAGNOSTICS ----
## cmdstanr doesn't automatically warn you when there are Rhat issues
fit2$cmdstan_diagnose()

# RANK PLOTS ----
# Find parameters that had poor mixing
summ <- fit2$draws() %>% summarise_draws()
poor_mixing <- summ %>%
  filter(rhat > 1.01 | ess_bulk < 400 | ess_tail < 400) %>%
  arrange(desc(rhat)) %>%
  pull(variable)

# Rank plots for poorly mixed parameters, if any
if (length(poor_mixing) == 0) {
  message("No parameters with poor mixing detected.")
} else {
  mcmc_rank_overlay(
    fit2$draws(),
    pars = poor_mixing
  )
}

# POSTERIOR RETRODICTIVE CHECKS ----
source("scripts/00_Function_plot_prc_hist.R")
plot_prc_hist(fit2, seedbank_df, max_count = 1500, n_bins = 100)
plot_prc_hist(fit2, seedbank_df, max_count = 100, n_bins = 100)
plot_prc_hist(fit2, seedbank_df, max_count = 20, n_bins = 10)
