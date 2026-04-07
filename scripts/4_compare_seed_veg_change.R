# Assess relationship between long-term above-ground vegetation changes and seed bank changes
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ---------------------------------------------------------------------
## Install cmdstanr/cdmstan if not installed already
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  #if cmdstanr is not installed
  install.packages(
    "cmdstanr",
    repos = c('https://stan-dev.r-universe.dev', getOption("repos"))
  ) #install cmdstanr
  cmdstanr::install_cmdstan() #install cmdstan
}
librarian::shelf(
  tidyverse,
  cmdstanr,
  tidybayes,
  bayesplot,
  posterior,
  ggdist,
  patchwork
)

# DATA -------------------------------------------------------------------------
## Download/load
source("scripts/source/download_data.R")
saddlegrid_habitat <- read.csv("data/saddptqd_xericmesic_categorization.csv") #categorization of long-term veg plots into mesic or xerix dry meadow

## Wrangle seed bank data
seedbank_df <- data_list$seedbank_composition.ac_hh.data %>%
  select(year:USDA_name, count) %>% #select relevant columns
  filter(
    substr(USDA_code, 1, 1) != 2 & USDA_code != "CAREX" & USDA_code != "POA"
  ) %>% #remove records not identified to species
  mutate(USDA_name = sub("\\s*var\\..*$", "", USDA_name)) %>%
  group_by(year, habitat, plot, USDA_code, USDA_name) %>%
  summarize(count = sum(count), .groups = "drop") %>%
  group_by(USDA_code, USDA_name) %>%
  mutate(count_scaled = scale(count)[, 1]) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_wider(
    names_from = year,
    values_from = count_scaled,
    names_prefix = "count"
  ) %>%
  mutate(change = count2023 - count1989)

## Wrangle long-term veg data
saddlegrid_df <- data_list$saddptqd.hh.data.csv %>%
  filter(
    year >= 1989 & year <= 2023, #same time frame as seed bank study
    plot %in% saddlegrid_habitat$plot,
    hit_type %in% c("top", "bottom") #middle and extra hits haven't been done across whole time period
  ) %>%
  mutate(USDA_name = sub("\\s*var\\..*$", "", USDA_name)) %>%
  group_by(year, plot, USDA_code, USDA_name) %>%
  summarize(count = length(hit_type), .groups = "drop") %>%
  complete(plot, USDA_code, year, fill = list(count = 0)) %>% #add records for species not seen in some years
  group_by(plot, USDA_code, USDA_name) %>%
  mutate(
    count_scaled = scale(count, center = FALSE)[, 1],
    count_scaled = ifelse(is.na(count_scaled), 0, count_scaled)
  ) %>%
  ungroup() %>%
  merge(saddlegrid_habitat, ., by = "plot") #add habitat information for later use

## Remove species-habitat combinations that weren't actually present in seed bank data
spp_to_include <- data_list$seedbank_composition.ac_hh.data %>%
  mutate(USDA_name = sub("\\s*var\\..*$", "", USDA_name)) %>%
  filter(
    substr(USDA_code, 1, 1) != 2 & USDA_code != "CAREX" & USDA_code != "POA"
  ) %>%
  filter(habitat == "mesic", count > 0) %>%
  distinct(USDA_name, habitat) %>%
  select(USDA_name, habitat) %>%
  rbind(
    data_list$seedbank_composition.ac_hh.data %>%
      mutate(USDA_name = sub("\\s*var\\..*$", "", USDA_name)) %>%
      filter(
        substr(USDA_code, 1, 1) != 2 & USDA_code != "CAREX" & USDA_code != "POA"
      ) %>%
      filter(habitat == "xeric", count > 0) %>%
      distinct(USDA_name, habitat) %>%
      select(USDA_name, habitat)
  )

### Filter seed bank data
seedf <- seedbank_df %>%
  filter(
    USDA_code %in% saddlegrid_df$USDA_code, #rm species that aren't in veg data
    paste0(USDA_name, habitat) %in%
      with(spp_to_include, paste0(USDA_name, habitat)) #remove species-habitat combos that werent present
  )

### Filter veg data to only include species now in seedf
vegf <- saddlegrid_df %>% filter(USDA_code %in% seedf$USDA_code)

## Create lookup data.frames for variables that will be scaled or turned into IDs for stan
year_lookup <- vegf %>%
  mutate(year_scaled = scale(year)[, 1]) %>%
  select(year, year_scaled) %>%
  distinct(year, year_scaled) %>%
  arrange(year)

USDA_lookup <- seedf %>%
  mutate(USDA_code_id = as.integer(factor(USDA_code))) %>%
  select(USDA_code_id, USDA_code, USDA_name) %>%
  distinct() %>%
  arrange(USDA_code_id)

veg_plot_lookup <- vegf %>%
  mutate(
    plot_id = as.integer(factor(plot)),
    habitat_id = habitat == "xeric"
  ) %>%
  select(plot_id, plot, habitat_id, habitat) %>%
  distinct() %>%
  arrange(plot_id)

seed_plot_lookup <- seedf %>%
  mutate(
    plot_id = as.integer(factor(plot)),
    habitat_id = habitat == "xeric"
  ) %>%
  select(plot_id, plot, habitat_id, habitat) %>%
  distinct() %>%
  arrange(plot_id)

# Put data into list for cmdstanr model
dlist <- list(
  # Vegetation data
  Nyear = nrow(vegf),
  veg_year = vegf$year - min(vegf$year),
  Kspecies = nrow(USDA_lookup),
  Nveg = nrow(vegf),
  veg_species = as.integer(factor(vegf$USDA_code)),
  veg_count = vegf$count_scaled,
  Nveg_plots = nrow(veg_plot_lookup),
  veg_plot = as.integer(factor(vegf$plot)),
  veg_habitat = veg_plot_lookup$habitat_id,

  # Seed data
  Nseed = nrow(seedf),
  seed_change = seedf$change / (2023 - 1989),
  seed_species = as.integer(factor(seedf$USDA_code)),
  Nseed_plots = nrow(seed_plot_lookup),
  seed_plot = as.integer(factor(seedf$plot)),
  seed_habitat = seed_plot_lookup$habitat_id
)

# MODEL ------------------------------------------------------------------------
options(mc.cores = ifelse(parallel::detectCores() > 4, 4, 2)) #set cores for parallel processing

## Compile
mod4 <- cmdstan_model("scripts/stan/compare_seed_veg_change.stan")

## Fit model
### will get warnings as model starts sampling at extreme values but fit is fine
fit4 <- mod4$sample(
  data = dlist,
  chains = 4,
  seed = 5336,
  adapt_delta = .99,
  max_treedepth = 15
)

## Check model diagnostics (model doesn't automatically warn you when there are Rhat issues)
fit4$cmdstan_diagnose()

## Posterior Predictive Checks
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
  "figures/figS5.pdf",
  ppc4,
  width = 10,
  height = 10,
  units = "in",
  dpi = 600
)

# FIGURES ----------------------------------------------------------------------
## Set plotting theme
theme_set(
  ggthemes::theme_tufte() +
    theme(
      panel.border = element_rect(fill = NA),
      axis.title.x = element_text(size = 8)
    )
)

## Get posterior draws for predicted seed change
seed_draws <- lapply(USDA_lookup$USDA_code_id, function(x) {
  tidy_draws(fit4) %>%
    mutate(
      USDA_code_id = x,
      mesic = get(paste0("a_seed[", x, "]")) +
        get(paste0("b_seed_mesic[", x, "]")) * # ← species-indexed
          get(paste0("b_veg[", x, ",1]")),
      xeric = get(paste0("a_seed[", x, "]")) +
        (get(paste0("b_seed_mesic[", x, "]")) + # ← species-indexed
          get(paste0("beta_seed_habitat_bveg[", x, "]"))) *
          get(paste0("b_veg[", x, ",2]")) +
        get(paste0("beta_seed_habitat[", x, "]")) * 1
    ) %>%
    select(.draw, USDA_code_id, mesic, xeric)
}) %>%
  do.call(rbind, .) %>%
  pivot_longer(
    cols = c(mesic, xeric),
    names_to = "habitat",
    values_to = "value"
  ) %>%
  mutate(USDA_name = USDA_lookup$USDA_name[USDA_code_id]) %>%
  select(.draw, habitat, USDA_name, value) %>%
  filter(
    paste0(USDA_name, habitat) %in%
      with(spp_to_include, paste0(USDA_name, habitat))
  )

## Get posterior draws for predicted veg change
veg_draws <- tidy_draws(fit4) %>%
  select(.draw, starts_with("b_veg")) %>%
  pivot_longer(cols = -.draw, names_to = "parameter", values_to = "value") %>%
  mutate(
    USDA_code_id = as.integer(str_extract(parameter, "(?<=\\[)\\d+")),
    USDA_name = USDA_lookup$USDA_name[USDA_code_id],
    habitat_id = as.integer(str_extract(parameter, "(?<=,)\\d+")),
    habitat = if_else(habitat_id == 1, "mesic", "xeric")
  ) %>%
  select(.draw, habitat, USDA_name, value) %>%
  filter(
    paste0(USDA_name, habitat) %in%
      with(spp_to_include, paste0(USDA_name, habitat))
  )

## Get posterior draws for relationship between veg and seed change
rel_draws <- lapply(USDA_lookup$USDA_code_id, function(x) {
  tidy_draws(fit4) %>%
    transmute(
      .draw,
      USDA_code_id = x,
      mesic = get(paste0("b_seed_mesic[", x, "]")),
      xeric = get(paste0("b_seed_mesic[", x, "]")) +
        get(paste0("beta_seed_habitat_bveg[", x, "]"))
    )
}) %>%
  do.call(rbind, .) %>%
  pivot_longer(
    cols = c(mesic, xeric),
    names_to = "habitat",
    values_to = "b_seed"
  ) %>%
  mutate(USDA_name = USDA_lookup$USDA_name[USDA_code_id]) %>%
  filter(
    paste0(USDA_name, habitat) %in%
      with(spp_to_include, paste0(USDA_name, habitat))
  )

## Figure out rank order for species on y-axis
spp_order <- seed_draws %>%
  filter(
    paste0(USDA_name, habitat) %in%
      with(spp_to_include, paste0(USDA_name, habitat))
  ) %>%
  group_by(USDA_name) %>%
  summarize(mean = mean(value), .groups = "drop") %>%
  arrange(mean) %>%
  pull(USDA_name) %>%
  as.character()

# Fig. 3A: Posteriors of predicted seed change
fig3A <- seed_draws %>%
  mutate(USDA_name = factor(as.character(USDA_name), levels = spp_order)) %>%
  ggplot(aes(x = value, y = fct_rev(USDA_name), fill = habitat)) +
  geom_hline(
    yintercept = spp_order,
    linetype = 2,
    color = "grey25",
    linewidth = 0.1
  ) +
  stat_slab(alpha = 0.5, normalize = "groups") +
  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  scale_x_continuous(
    name = "Scaled Change per Year in Seed Bank",
    limits = c(-0.2, 0.2)
  ) +
  scale_y_discrete(name = "Species") +
  scale_fill_manual(
    values = c("mesic" = "green4", "xeric" = "tan4"),
    name = "Habitat"
  ) +
  theme(axis.text.y = element_text(face = "italic"))

# Fig. 3A: Posteriors of predicted veg change
fig3B <- veg_draws %>%
  mutate(USDA_name = factor(as.character(USDA_name), levels = spp_order)) %>%
  ggplot(aes(x = value, y = fct_rev(USDA_name), fill = habitat)) +
  geom_hline(
    yintercept = spp_order,
    linetype = 2,
    color = "grey25",
    linewidth = 0.1
  ) +
  stat_slab(alpha = 0.5, normalize = "groups") +
  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  scale_x_continuous(name = "Scaled Change per Year in Vegetation") +
  scale_y_discrete(name = "Species") +
  scale_fill_manual(
    values = c("mesic" = "green4", "xeric" = "tan4"),
    name = "Habitat"
  ) +
  theme(axis.text.y = element_text(face = "italic"))

# Fig. 3A: Posteriors of relationship between seed and veg change
fig3C <- rel_draws %>%
  ggplot(aes(x = b_seed, fill = habitat)) +
  stat_slab(alpha = 0.5, normalize = "groups") +
  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  scale_x_continuous(
    name = "Relationship between Seed Change and Veg Change",
    limits = c(-20, 20)
  ) +
  scale_fill_manual(
    values = c("mesic" = "green4", "xeric" = "tan4"),
    guide = "none"
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

fig3 <- ((fig3A + fig3B + plot_layout(axes = "collect")) / fig3C) +
  plot_layout(guides = "collect")

ggsave(
  "figures/fig3.pdf",
  fig3,
  width = 7.5,
  height = 5,
  units = "in",
  dpi = 600
)
ggsave(
  "figures/fig3.png",
  fig3,
  width = 10,
  height = 4,
  units = "in",
  dpi = 600,
  bg = "white"
)


# MODEL w/o HABITAT ------------------------------------------------------------
## Compile
mod4_wo_habitat <- cmdstan_model(
  "scripts/stan/compare_seed_veg_change_wo_habitat.stan"
)

## Fit
### will get warnings as model starts sampling at extreme values but fit is fine
fit4_wo_habitat <- mod4_wo_habitat$sample(
  data = dlist,
  chains = 4,
  seed = 5336,
  adapt_delta = .95
)

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

## FIGURES ----------------------------------------------------------------------
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
