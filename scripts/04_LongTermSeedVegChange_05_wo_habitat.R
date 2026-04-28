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
  "scripts/stan_models/compare_seed_veg_change_wo_habitat.stan"
)

# Fit
## will get warnings as model starts sampling at extreme values but fit is fine
fit4_wo_habitat <- mod4_wo_habitat$sample(
  data = dlist,
  chains = 4,
  parallel_chains = ifelse(parallel::detectCores() > 4, 4, 2),
  seed = 5336,
  adapt_delta = 0.95
)

# CHECKS ----
# Check model diagnostics (model doesn't always automatically warn you when there are issues)
fit4_wo_habitat$cmdstan_diagnose()

## Posterior Predictive Checks
### Extract draws
seed_rep_wo_habitat <- fit4_wo_habitat$draws(
  "seed_change_rep",
  format = "matrix"
)
veg_rep_wo_habitat <- fit4_wo_habitat$draws("veg_count_rep", format = "matrix")
y_seed <- dlist$seed_change
y_veg <- dlist$veg_count

### Grouped by habitat and species
ppc_stat_grouped(
  y_seed,
  seed_rep_wo_habitat,
  group = USDA_lookup$USDA_code[dlist$seed_species],
  stat = "mean"
)
ppc_stat_grouped(
  y_seed,
  seed_rep_wo_habitat,
  group = USDA_lookup$USDA_code[dlist$seed_species],
  stat = "sd"
)
ppc_stat_grouped(
  y_veg,
  veg_rep_wo_habitat,
  group = USDA_lookup$USDA_code[dlist$veg_species],
  stat = "mean"
)
ppc_stat_grouped(
  y_veg,
  veg_rep_wo_habitat,
  group = USDA_lookup$USDA_code[dlist$veg_species],
  stat = "sd"
)

# VIZ ----
theme_set(
  ggthemes::theme_tufte() + theme(panel.border = element_rect(fill = NA))
) #set theme for plotting

### Get posterior draws for predicted seed change
seed_draws_wo_habitat <- lapply(USDA_lookup$USDA_code_id, function(x) {
  tidy_draws(fit4_wo_habitat) %>%
    mutate(
      USDA_code_id = x,
      epred = get(paste0("a_seed[", x, "]")) +
        get("b_seed") * get(paste0("b_veg[", x, "]"))
    ) %>%
    select(USDA_code_id, epred)
}) %>%
  do.call(rbind, .) %>%
  mutate(
    USDA_name = factor(USDA_lookup$USDA_name[USDA_code_id], levels = spp_order)
  )

### Get posterior draws for predicted veg change
veg_draws_wo_habitat <- tidy_draws(fit4_wo_habitat) %>%
  select(starts_with("b_veg")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "b_veg"
  ) %>%
  mutate(
    USDA_code_id = as.integer(str_extract(parameter, "(?<=\\[)\\d+")),
    USDA_name = factor(USDA_lookup$USDA_name[USDA_code_id], levels = spp_order)
  )

### Get posterior draws for relationship between veg and seed change
rel_draws_wo_habitat <- tidy_draws(fit4_wo_habitat) %>%
  select(b_seed)

# Posterior distributions for change in seed abundance for each species
figS3A <- seed_draws_wo_habitat %>%
  ggplot(aes(x = epred, y = fct_rev(USDA_name))) +
  geom_hline(
    yintercept = spp_order,
    linetype = 2,
    color = "grey25",
    linewidth = 0.1
  ) +
  stat_slab(fill = "black", normalize = "groups") +
  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  scale_x_continuous(
    name = "Scaled Change per Year in Seed Bank",
    limits = c(-0.2, 0.2)
  ) +
  scale_y_discrete(name = "Species") +
  theme(axis.text.y = element_text(face = "italic"))

# Posterior distributions for change in veg cover for each species
figS3B <- veg_draws_wo_habitat %>%
  ggplot(aes(x = b_veg, y = fct_rev(USDA_name))) +
  geom_hline(
    yintercept = spp_order,
    linetype = 2,
    color = "grey25",
    linewidth = 0.1
  ) +
  stat_slab(fill = "black", normalize = "groups") +
  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  scale_x_continuous(name = "Scaled Change per Year in Vegetation") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Posterior distribution for relationship between seed and veg change
figS3C <- rel_draws_wo_habitat %>%
  ggplot(aes(x = b_seed)) +
  stat_slab(fill = "black", normalize = "groups") +
  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  scale_x_continuous(name = "Relationship between Seed Change and Veg Change") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

figS3 <- (figS3A | figS3B + plot_layout(axes = "collect")) /
  figS3C +
  plot_layout(guides = "collect") +
  plot_annotation(tag_level = "A")

# Write Figure 3
ggsave(
  "figures/figS3.pdf",
  figS3,
  width = 7.5,
  height = 5,
  units = "in",
  dpi = 600
)
