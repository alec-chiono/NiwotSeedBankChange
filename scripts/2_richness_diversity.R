# Estimate richness, evenness, and Hill's N1 and change in these indices over time for soil community
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ---------------------------------------------------------------------
## Install cmdstanr/cdmstan if not installed already
if (!requireNamespace("cmdstanr", quietly=TRUE)) { #if cmdstanr is not installed
  install.packages("cmdstanr", repos=c('https://stan-dev.r-universe.dev', getOption("repos"))) #install cmdstanr
  cmdstanr::install_cmdstan() #install cmdstan
}
## Load packages
librarian::shelf(tidyverse, cmdstanr, tidybayes, bayesplot, posterior, loo, ggdist, patchwork)

# DATA -------------------------------------------------------------------------
## Download
source("scripts/source/download_data.R")

## Wrangle
seedbank_df <- data_list$seedbank_composition.ac_hh.data %>%
  select(year:USDA_code, count) %>%  #select relevant columns
  filter(substr(USDA_code,1, 1)!=2 & USDA_code!="CAREX" & USDA_code!="POA") %>%  #remove records not identified to species
  group_by(year, habitat, plot, USDA_code) %>%
  summarize(count=sum(count), .groups="drop")

## Put data into list for cmdstanr model
dlist <- list(
  N_obs=nrow(seedbank_df),
  N_species=length(unique(seedbank_df$USDA_code)),
  N_plots=length(unique(seedbank_df$plot)),
  N_years=length(unique(seedbank_df$year)),
  N_habitats=length(unique(seedbank_df$habitat)),
  count=seedbank_df$count,
  species=as.integer(factor(seedbank_df$USDA_code)),
  plot=as.integer(factor(seedbank_df$plot)),
  year_idx=as.integer(factor(seedbank_df$year)),
  habitat_idx=as.integer(factor(seedbank_df$habitat))
)

# MODEL ------------------------------------------------------------------------
options(mc.cores=ifelse(parallel::detectCores()>4, 4, 2)) #set cores for parallel processing

## Compile model
mod2 <- cmdstan_model("scripts/stan/richness_diversity.stan")

## Fit model
### will get warnings as model starts sampling at extreme values but fit is fine
fit2 <- mod2$sample(
  data=dlist,
  chains=4,
  seed=4879523,
  iter_sampling=1000,
  adapt_delta=0.99
)

## Check model diagnostics (model doesn't always automatically warn you when there are issues)
fit2$cmdstan_diagnose()

## Posterior Predictive Check
log_lam <- fit2$draws("log_lambda_gq", format="matrix")
phi_draws <- fit2$draws("phi", format="matrix")
yrep <- matrix(
  rnbinom(
    nrow(log_lam)*dlist$N_obs,
    mu=exp(as.vector(log_lam)),
    size=as.vector(phi_draws)),
  nrow=nrow(log_lam)
  )
ppc_stat(dlist$count, yrep, stat = \(y) mean(y == 0)) #proportion of zeroes
ppc_stat(dlist$count, yrep, stat=\(y)var(y)/mean(y)) #dispersion ratio
ppc_stat(dlist$count, yrep_samp, stat = \(y) quantile(y, 0.9)) #90th percentile (i.e. tail-heaviness)
ppc_stat(dlist$count, yrep_samp, stat = "median") #mean
ppc_stat(dlist$count, yrep, stat="sd") #std dev

lv <- fit2$loo(save_psis = TRUE)
ppc_loo_pit_qq(dlist$count, yrep, psis_object = lv$psis_object)
plot(lv)  # Visual of k values by observation
seedbank_df[which(lv$diagnostics$pareto_k > 0.7), ]  # Which observation have high k


# FIGURES ----------------------------------------------------------------------
## Set plotting theme
theme_set(ggthemes::theme_tufte() + theme(panel.border=element_rect(fill=NA)))

## Retrieve posteriors for seed counts
post2_counts <-
  tidy_draws(fit2) %>%
  select(.draw, starts_with("predicted_seed_counts")) %>%
  pivot_longer(-.draw,
               names_to = c("metric", "year", "habitat", "species"),
               names_pattern = "(predicted_seed_counts)\\[(\\d+),(\\d+),(\\d+)\\]",
               values_to = "raw_count") %>%
  group_by(species) %>%
  mutate(scaled_count = raw_count / sd(raw_count)) %>%
  group_by(.draw, year, habitat) %>%
  summarize(scaled_total = sum(scaled_count), .groups="drop") %>%
  mutate(
    metric = "Scaled Total Seeds",
    year = factor(year, levels = 1:dlist$N_years, labels = sort(unique(seedbank_df$year))),
    habitat = factor(habitat, levels = 1:dlist$N_habitats, labels = sort(unique(seedbank_df$habitat))),
    value = scaled_total
  ) %>%
  select(.draw, metric, year, habitat, value)

## Retrieve posteriors for richness, evenness, and diversity indices
post2_indices <-
  tidy_draws(fit2) %>%
  select(.draw, starts_with("richness"), starts_with("evenness"), starts_with("hill_N1")) %>%
  pivot_longer(-.draw,
               names_to=c("metric", "year", "habitat"),
               names_pattern="(predicted_seed_counts|richness|evenness|hill_N1)\\[(\\d+),(\\d+)]",
               values_to="value") %>%
  mutate(
    metric=factor(metric,
                  levels=c("richness", "evenness", "hill_N1"),
                  labels=c("Species Richness", "Pielou's Evenness", "Hill's Diversity Index")),
    year=factor(year, levels=1:dlist$N_years, labels=sort(unique(seedbank_df$year))),
    habitat=factor(habitat, levels=1:dlist$N_habitats, labels=sort(unique(seedbank_df$habitat)))
  )

## Collate posteriors together
post2 <- bind_rows(post2_indices, post2_counts) %>%
  mutate(metric=factor(metric, levels=c("Total Seeds Scaled", "Species Richness", "Pielou's Evenness", "Hill's Diversity Index")))

## Fig. 1B: Estimates for each year and habitat
fig1B <- post2 %>%
  mutate(
    y_dots = if_else(metric == "Species Richness", value, NA),
    y_slab = if_else(metric != "Species Richness", value, NA)
  ) %>%
  ggplot(aes(x=year)) +
  stat_dots(aes(y=y_dots), fill="NA",color="black", linewidth=0.1, normalize="panels") +
  stat_slab(aes(y=y_slab), fill="black",color="black", linewidth=0.1, normalize="panels") +
  facet_grid(metric ~ habitat, scales="free_y") +
  scale_y_continuous(name="Estimated Value")

## Fig. 1C Estimates for change in indices over time
fig1C <- post2 %>%
  arrange(.draw, metric, habitat, year) %>%
  group_by(.draw, metric, habitat) %>%
  mutate(diff=lead(value)-value, .groups="drop") %>%
  filter(!is.na(diff)) %>%
  mutate(
    y_dots = if_else(metric == "Species Richness", diff, NA),
    y_slab = if_else(metric != "Species Richness", diff, NA)
  ) %>%
  ggplot() +
  stat_dots(aes(y=y_dots), fill="NA",color="black", linewidth=0.1, normalize="panels") +
  stat_slab(aes(y=y_slab), fill="black",color="black", linewidth=0.1, normalize="panels") +
  geom_hline(yintercept=0, linetype="dashed", color="red") +
  facet_grid(metric ~ habitat, scales="free_y") +
  scale_y_continuous(name="Difference between Years") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

## Get ggplot obj for Fig. 1A
source("scripts/1_seed_community_RDA.R")

## Collate Figure 1
fig1 <- fig1A / (fig1B + fig1C) + plot_annotation(tag_levels = 'A')

## Write Figure 1
ggsave("figures/fig1.pdf", fig1, width=7.5, height=10, units="in", dpi=600)
