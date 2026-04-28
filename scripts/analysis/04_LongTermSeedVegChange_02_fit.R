# Fit Bayesian model evaluating relationship between long-term above-ground vegetation changes and seed bank changes
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ----
# Install cmdstanr/cdmstan if not installed already
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  #if cmdstanr is not installed
  install.packages(
    "cmdstanr",
    repos = c('https://stan-dev.r-universe.dev', getOption("repos"))
  ) #install cmdstanr
  cmdstanr::install_cmdstan() #install cmdstan
}
library(cmdstanr)

# DATA ----
source("scripts/analysis/04_LongTermSeedVegChange_01_prep.R")

# MODEL ----
# Compile
mod4 <- cmdstan_model("scripts/stan_models/compare_seed_veg_change.stan")

# Fit
## may get warnings as model starts sampling at extreme values but fit is fine
fit4 <- mod4$sample(
  data = dlist,
  chains = 4,
  parallel_chains = ifelse(parallel::detectCores() > 4, 4, 2),
  seed = 5336,
  adapt_delta = 0.99,
  max_treedepth = 15
)
