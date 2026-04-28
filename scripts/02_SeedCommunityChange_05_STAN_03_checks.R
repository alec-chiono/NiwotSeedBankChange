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

# PPC ----
## Posterior Predictive Checka
log_lam <- fit2$draws("log_lambda_gq", format = "matrix")
phi_draws <- fit2$draws("phi", format = "matrix")
yrep <- matrix(
  rnbinom(
    nrow(log_lam) * dlist$N_obs,
    mu = exp(as.vector(log_lam)),
    size = as.vector(phi_draws)
  ),
  nrow = nrow(log_lam)
)

ppc2A <- ppc_stat(dlist$count, yrep, stat = \(y) mean(y == 0)) #proportion of zeroes
ppc2B <- ppc_stat(dlist$count, yrep, stat = \(y) var(y) / mean(y)) #dispersion ratio
ppc2C <- ppc_stat(dlist$count, yrep, stat = \(y) quantile(y, 0.9)) #90th percentile (i.e. tail-heaviness)
ppc2D <- ppc_stat(dlist$count, yrep, stat = "median") #mean
ppc2E <- ppc_stat(dlist$count, yrep, stat = "sd") #std dev

ppc2 <- ppc2A /
  ppc2B /
  ppc2C /
  ppc2D /
  ppc2E +
  plot_annotation(tag_levels = "A")

ggsave(
  "figures/figS4.pdf",
  ppc2,
  width = 10,
  height = 10,
  units = "in",
  dpi = 600
)

lv <- fit2$loo(save_psis = TRUE)
ppc_loo_pit_qq(dlist$count, yrep, psis_object = lv$psis_object)
plot(lv) # Visual of k values by observation
seedbank_df[which(lv$diagnostics$pareto_k > 0.7), ] # Which observation have high k
