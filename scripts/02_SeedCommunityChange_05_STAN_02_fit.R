# Fit Bayesian model for estimating richness, evenness, and diversity
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
source("scripts/02_SeedCommunityChange_05_STAN_01_prep.R")

# MODEL ----
# Compile
mod2 <- cmdstan_model("scripts/stan_models/community_indices.stan")

# Fit
## may get warnings as model starts sampling at extreme values but fit is fine
fit2 <- mod2$sample(
  data = dlist,
  chains = 4,
  parallel_chains = ifelse(parallel::detectCores() > 4, 4, 2),
  seed = 11001001
)
