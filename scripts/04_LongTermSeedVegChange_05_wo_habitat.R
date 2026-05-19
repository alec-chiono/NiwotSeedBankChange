# Replicate Bayesian model evaluating relationship between long-term above-ground vegetation changes and seed bank changes but ignoring habitat
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
library(bayesplot)
library(posterior)
library(patchwork)
library(tidybayes)
library(ggdist)

# DATA ----
source("scripts/04_LongTermSeedVegChange_01_prep.R")

# MODEL ----
# Compile
mod4_wo_habitat <- cmdstan_model(
  "scripts/stan_models/longterm_change_wo_habitat.stan"
)

# Fit
## will get warnings as model starts sampling at extreme values but fit is fine
fit4_wo_habitat <- mod4_wo_habitat$sample(
  data = dlist,
  chains = 4,
  parallel_chains = ifelse(parallel::detectCores() > 4, 4, 2),
  seed = 5336,
  adapt_delta = 0.99
)

# DIAGNOSTICS ----
## cmdstanr doesn't automatically warn you when there are Rhat issues
fit4_wo_habitat$cmdstan_diagnose()

# RANK PLOTS ----
# Find parameters that had poor mixing
summ <- fit4_wo_habitat$draws() %>% summarise_draws()
poor_mixing <- summ %>%
  filter(rhat > 1.01 | ess_bulk < 400 | ess_tail < 400) %>%
  arrange(desc(rhat)) %>%
  pull(variable)

# Rank plots for poorly mixed parameters, if any
if (length(poor_mixing) == 0) {
  message("No parameters with poor mixing detected.")
} else {
  mcmc_rank_overlay(
    fit4_wo_habitat$draws(),
    pars = poor_mixing
  )
}

# VIZ ----
theme_set(
  ggthemes::theme_tufte() + theme(panel.border = element_rect(fill = NA))
) #set theme for plotting

source("scripts/04_LongTermSeedVegChange_04_viz.R") #get spp_order

# Get posterior draws for predicted seed change
seed_draws_wo_habitat <- lapply(USDA_lookup$USDA_code_id, function(x) {
  tidy_draws(fit4_wo_habitat) %>%
    mutate(
      USDA_code_id = x,
      epred = get(paste0("a_seed[", x, "]")) +
        get(paste0("b_seed[", x, "]")) * # species-specific slope; was get("b_seed")
          get(paste0("b_veg[", x, "]"))
    ) %>%
    select(USDA_code_id, epred)
}) %>%
  do.call(rbind, .) %>%
  mutate(
    Species = factor(USDA_lookup$Species[USDA_code_id], levels = spp_order)
  )

# Get posterior draws for predicted veg change
veg_draws_wo_habitat <- tidy_draws(fit4_wo_habitat) %>%
  select(matches("^b_veg\\[\\d+\\]$")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "b_veg"
  ) %>%
  mutate(
    USDA_code_id = as.integer(str_extract(parameter, "(?<=\\[)\\d+")),
    Species = factor(USDA_lookup$Species[USDA_code_id], levels = spp_order)
  )

# Get posterior draws for the community-level vegetation–seed linkage
rel_draws_wo_habitat <- tidy_draws(fit4_wo_habitat) %>%
  select(mu_b_seed) # was: select(b_seed)

# Posterior distributions for change in seed abundance for each species
figS3A <- seed_draws_wo_habitat %>%
  ggplot(aes(x = epred, y = fct_rev(Species))) +
  geom_hline(
    yintercept = spp_order,
    linetype = 2,
    color = "grey25",
    linewidth = 0.1
  ) +
  stat_slab(fill = "black", normalize = "groups") +
  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  scale_x_continuous(
    name = "Annual change in seed bank abundance",
    limits = c(-0.2, 0.2)
  ) +
  scale_y_discrete(name = "Species") +
  theme(axis.text.y = element_text(face = "italic"))

# Posterior distributions for change in veg cover for each species
figS3B <- veg_draws_wo_habitat %>%
  ggplot(aes(x = b_veg, y = fct_rev(Species))) +
  geom_hline(
    yintercept = spp_order,
    linetype = 2,
    color = "grey25",
    linewidth = 0.1
  ) +
  stat_slab(fill = "black", normalize = "groups") +
  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  scale_x_continuous(
    name = "Annual change in vegetation abundance",
    limits = c(-0.05, 0.05)
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Posterior distribution for the community-level vegetation–seed linkage
figS3C <- rel_draws_wo_habitat %>%
  ggplot(aes(x = mu_b_seed)) + # was: aes(x = b_seed)
  stat_slab(fill = "black", normalize = "groups") +
  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  scale_x_continuous(
    name = "Relationship between seed bank change and vegetation change"
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

figS3 <- (figS3A | figS3B + plot_layout(axes = "collect")) /
  figS3C +
  plot_layout(guides = "collect") +
  plot_annotation(tag_level = "A")

# Write Figure S3
ggsave(
  "figures/FigureS3.pdf",
  figS3,
  width = 7.5,
  height = 5,
  units = "in",
  dpi = 600
)
