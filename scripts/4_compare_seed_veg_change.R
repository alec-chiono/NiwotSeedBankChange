# Assess relationship between long-term above-ground vegetation changes and seed bank changes
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ---------------------------------------------------------------------
## Install cmdstanr/cdmstan if not installed already
if (!requireNamespace("cmdstanr", quietly=TRUE)) { #if cmdstanr is not installed
  install.packages("cmdstanr", repos=c('https://stan-dev.r-universe.dev', getOption("repos"))) #install cmdstanr
  cmdstanr::install_cmdstan() #install cmdstan
}
librarian::shelf(tidyverse, cmdstanr, tidybayes, bayesplot, ggdist, patchwork)

# DATA -------------------------------------------------------------------------
## Download/load
source("Zenodo_archiving/scripts/source/download_data.R")
saddlegrid_habitat <- read.csv("Zenodo_archiving/data/saddptqd_xericmesic_categorization.csv") #categorization of long-term veg plots into mesic or xerix dry meadow

## Wrangle seed bank data
seedbank_df <- data_list$seedbank_composition.ac_hh.data  %>%
  select(year:USDA_name, count) %>%  #select relevant columns
  filter(substr(USDA_code,1, 1)!=2 & USDA_code!="CAREX" & USDA_code!="POA") %>%  #remove records not identified to species
  group_by(year, habitat, plot, USDA_code, USDA_name) %>%
  summarize(count=sum(count), .groups="drop") %>%
  group_by(USDA_code, USDA_name) %>%
  mutate(count_scaled=scale(count)[,1]) %>% ungroup() %>%
  select(-count) %>%
  pivot_wider(names_from=year, values_from=count_scaled, names_prefix="count") %>%
  mutate(change=count2023-count1989)

## Wrangle long-term veg data
saddlegrid_df <- data_list$saddptqd.hh.data.csv %>%
  filter(
    plot%in%saddlegrid_habitat$plot,
    hit_type%in%c("top","bottom") #middle and extra hits haven't been done across whole time period
  ) %>%
  group_by(year, plot, USDA_code) %>%
  summarize(count=length(hit_type), .groups="drop") %>%
  complete(plot, USDA_code, year, fill=list(count=0)) %>% #add records for species not seen in some years
  group_by(plot, USDA_code) %>%
  mutate(count_scaled=scale(count, center=FALSE)[,1],
         count_scaled=ifelse(is.na(count_scaled), 0, count_scaled)) %>% ungroup() %>%
  merge(saddlegrid_habitat, ., by="plot") #add habitat information for later use

## Remove species-habitat combinations that weren't actually present in seed bank data
### Figure out species present in just mesic seed bank
spp_to_include <- data_list$seedbank_composition.ac_hh.data   %>%
  filter(substr(USDA_code,1, 1)!=2 & USDA_code!="CAREX" & USDA_code!="POA") %>%
  filter(habitat == "mesic", count>0) %>%
  distinct(USDA_name, habitat) %>%
  select(USDA_name, habitat) %>%
  rbind(data_list$seedbank_composition.ac_hh.data   %>%
          filter(substr(USDA_code,1, 1)!=2 & USDA_code!="CAREX" & USDA_code!="POA") %>%
          filter(habitat == "xeric", count>0) %>%
          distinct(USDA_name, habitat) %>%
          select(USDA_name, habitat)
        )

### Filter seed bank data
seedf <- seedbank_df %>%
  filter(paste0(USDA_name, habitat) %in% with(spp_to_include, paste0(USDA_name, habitat)))

### Filter veg data to only include speceis in seed bank
vegf <- saddlegrid_df %>% filter(USDA_code%in%seedf$USDA_code)


## Create lookup data.frames for variables that will be scaled or turned into IDs for stan
year_lookup <- vegf %>%
  mutate(year_scaled=scale(year)[,1]) %>%
  select(year, year_scaled) %>%
  distinct(year, year_scaled) %>%
  arrange(year)

USDA_lookup <- seedf %>%
  mutate(USDA_code_id=as.integer(factor(USDA_code))) %>%
  select(USDA_code_id, USDA_code, USDA_name) %>%
  distinct() %>%
  arrange(USDA_code_id)

veg_plot_lookup <- vegf %>%
  mutate(
    plot_id=as.integer(factor(plot)),
    habitat_id=habitat=="xeric"
  ) %>%
  select(plot_id, plot, habitat_id, habitat) %>%
  distinct() %>%
  arrange(plot_id)

seed_plot_lookup <- seedf %>%
  mutate(
    plot_id=as.integer(factor(plot)),
    habitat_id=habitat=="xeric"
  ) %>%
  select(plot_id, plot, habitat_id, habitat) %>%
  distinct() %>%
  arrange(plot_id)

# Put data into list for cmdstanr model
dlist <- list(
  # Vegetation data
  Nyear=nrow(vegf),
  veg_year=vegf$year-min(vegf$year),
  Kspecies=nrow(USDA_lookup),
  Nveg=nrow(vegf),
  veg_species=as.integer(factor(vegf$USDA_code)),
  veg_count=vegf$count_scaled,
  Nveg_plots=nrow(veg_plot_lookup),
  veg_plot=as.integer(factor(vegf$plot)),
  veg_habitat=veg_plot_lookup$habitat_id,

  # Seed data
  Nseed=nrow(seedf),
  seed_change=seedf$change/(2023-1989),
  seed_species=as.integer(factor(seedf$USDA_code)),
  Nseed_plots=nrow(seed_plot_lookup),
  seed_plot=as.integer(factor(seedf$plot)),
  seed_habitat=seed_plot_lookup$habitat_id
)

# MODEL ------------------------------------------------------------------------
options(mc.cores=ifelse(parallel::detectCores()>4, 4, 2)) #set cores for parallel processing

## Compile
mod4 <- cmdstan_model("Zenodo_archiving/scripts/stan/compare_seed_veg_change_w_habitat.stan")

## Fit model
### will get warnings as model starts sampling at extreme values but fit is fine
fit4 <- mod4$sample(
  data=dlist,
  chains=4,
  seed=5336,
  adapt_delta=.99,
  max_treedepth=15
)

## Check model diagnostics (model doesn't always automatically warn you when there are issues)
fit4$cmdstan_diagnose()

## Posterior Predictive Checks

# FIGURES ----------------------------------------------------------------------
## Set plotting theme
theme_set(ggthemes::theme_tufte() + theme(panel.border=element_rect(fill=NA)))

## Get posterior draws for predicted seed change
seed_draws <- lapply(USDA_lookup$USDA_code_id, function(x) {
  tidy_draws(fit4) %>%
    mutate(
      USDA_code_id=x,
      # Expected seed change for MESIC habitat (habitat=0 -> b_veg col 1)
      mesic=get(paste0("a_seed[", x, "]")) +
        b_seed_mesic * get(paste0("b_veg[", x, ",1]")) +
        get(paste0("beta_seed_habitat[", x, "]")) * 0 +
        get(paste0("beta_seed_habitat_bveg[", x, "]")) * 0 * get(paste0("b_veg[", x, ",1]")),
      # Expected seed change for XERIC habitat (habitat=1 -> b_veg col 2)
      xeric=get(paste0("a_seed[", x, "]")) +
        (b_seed_mesic + get(paste0("beta_seed_habitat_bveg[", x, "]"))) * get(paste0("b_veg[", x, ",2]")) +
        get(paste0("beta_seed_habitat[", x, "]")) * 1
    ) %>%
    select(USDA_code_id, mesic, xeric)
}) %>%
  do.call(rbind, .) %>%
  pivot_longer(cols=c(mesic, xeric), names_to="habitat", values_to="value") %>%
  mutate(USDA_name=USDA_lookup$USDA_name[USDA_code_id]) %>%
  select(habitat, USDA_name, value) %>%
  filter(paste0(USDA_name, habitat) %in% with(spp_to_include, paste0(USDA_name, habitat)))

## Get posterior draws for predicted veg change
veg_draws <- fit4$draws("b_veg", format="df") %>%
  select(starts_with("b_veg")) %>%
  pivot_longer(cols=everything(), names_to="parameter", values_to="value") %>%
  mutate(
    USDA_code_id=as.integer(str_extract(parameter, "(?<=\\[)\\d+")),
    USDA_name=USDA_lookup$USDA_name[USDA_code_id],
    habitat_id=as.integer(str_extract(parameter, "(?<=,)\\d+")),
    habitat=if_else(habitat_id == 1, "mesic", "xeric")
  ) %>%
  select(habitat, USDA_name, value) %>%
  filter(paste0(USDA_name, habitat) %in% with(spp_to_include, paste0(USDA_name, habitat)))

## Get posterior draws for relationship between veg and seed change
rel_draws <- tidy_draws(fit4) %>%
  mutate(b_seed_xeric=b_seed_mesic + rowMeans(across(starts_with("beta_seed_habitat_bveg")))) %>%
  select(b_seed_mesic, b_seed_xeric) %>%
  pivot_longer(
    cols=everything(),
    names_to="habitat",
    values_to="b_seed"
  ) %>%
  mutate(habitat=recode(habitat, b_seed_mesic ="mesic", b_seed_xeric ="xeric"))

## Figure out rank order for species on y-axis
spp_order <- seed_draws %>%
  filter(paste0(USDA_name, habitat) %in% with(spp_to_include, paste0(USDA_name, habitat))) %>%
  group_by(USDA_name) %>%
  summarize(mean=mean(value), .groups="drop") %>%
  arrange(mean) %>%
  pull(USDA_name) %>% as.character()

# Fig. 3A: Posteriors of predicted seed change
fig3A <- seed_draws %>%
  mutate(USDA_name=factor(as.character(USDA_name), levels=spp_order)) %>%
  ggplot(aes(x=value, y=(USDA_name), fill=habitat)) +
  geom_hline(yintercept=spp_order, linetype=2, color="grey25", linewidth=0.1) +
  stat_slab(alpha=0.75, normalize="groups") +
  geom_vline(xintercept=0, linetype=2, color="red") +
  scale_x_continuous(name="Change") +
  scale_y_discrete(name="Species") +
  scale_fill_manual(values=c("mesic"="green", "xeric"="tan3"), name="Habitat") +
  theme(axis.text.y=element_text(face="italic"))

# Fig. 3A: Posteriors of predicted veg change
fig3B <- veg_draws %>%
  mutate(USDA_name=factor(as.character(USDA_name), levels=spp_order)) %>%
  ggplot(aes(x=value, y=(USDA_name), fill=habitat)) +
  geom_hline(yintercept=spp_order, linetype=2, color="grey25", linewidth=0.1) +
  stat_slab(alpha=0.75, normalize="groups") +
  geom_vline(xintercept=0, linetype=2, color="red") +
  scale_x_continuous(name="Change") +
  scale_y_discrete(name="Species") +
  scale_fill_manual(values=c("mesic"="green", "xeric"="tan3"), name="Habitat") +
  theme(axis.text.y=element_text(face="italic"))

# Fig. 3A: Posteriors of relationship between seed and veg change
fig3C <- rel_draws %>%
  ggplot(aes(x=b_seed, fill=habitat)) +
  stat_slab(alpha=0.6, normalize="groups") +
  geom_vline(xintercept=0, linetype=2, color="red") +
  scale_x_continuous(name="Seed Change ~ Veg Change") +
  scale_fill_manual(values=c("mesic"="green", "xeric"="tan3"), guide="none") +
  theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
  )

fig3 <- ((fig3A + fig3B + plot_layout(axes="collect")) / fig3C) + plot_layout(guides="collect")

ggsave("Zenodo_archiving/figures/fig3.pdf", fig3, width=7.5, height=5, units="in", dpi=600)


# Model without habitat --------------------------------------------------------
## Wrangle data
## Seed bank data
seedbank_df <- data_list$seedbank_composition.ac_hh.data  %>%
  select(year:USDA_name, count) %>%  #select relevant columns
  filter(substr(USDA_code,1, 1)!=2 & USDA_code!="CAREX" & USDA_code!="POA") %>%  #remove records not identified to species
  group_by(year, habitat, plot, USDA_code, USDA_name) %>%
  summarize(count=sum(count), .groups="drop") %>%
  group_by(USDA_code, USDA_name) %>%
  mutate(count_scaled=scale(count)[,1]) %>% ungroup() %>%
  select(-count) %>%
  pivot_wider(names_from=year, values_from=count_scaled, names_prefix="count") %>%
  mutate(change=count2023-count1989)

## Long-term veg data
saddlegrid_df <- data_list$saddptqd.hh.data.csv %>%
  filter(
    plot%in%saddlegrid_habitat$plot,
    hit_type%in%c("top","bottom") #middle and extra hits haven't been done across whole time period
  ) %>%
  group_by(year, plot, USDA_code) %>%
  summarize(count=length(hit_type), .groups="drop") %>%
  complete(plot, USDA_code, year, fill=list(count=0)) %>% #add records for species not seen in some years
  group_by(plot, USDA_code) %>%
  mutate(count_scaled=scale(count, center=FALSE)[,1],
         count_scaled=ifelse(is.na(count_scaled), 0, count_scaled)) %>% ungroup() %>%
  merge(saddlegrid_habitat, ., by="plot") #add habitat information for later use

## Pare down both data sets to include only species present in both data sets
veg_USDAs <- saddlegrid_df %>% distinct(USDA_code)
seed_USDAs <- seedbank_df %>% distinct(USDA_code)

### Species present in BOTH data.frames
both_USDA <- inner_join(veg_USDAs, seed_USDAs, by=c("USDA_code"))$USDA_code

## Filter both data.frames to include only species present in both
vegf <- saddlegrid_df %>% filter(USDA_code%in%both_USDA)
seedf <- seedbank_df %>% filter(USDA_code%in%both_USDA) %>% arrange(USDA_code)

## Create lookup data.frames for variables that will be scaled or turned into IDs for stan
year_lookup <- vegf %>%
  mutate(year_scaled=scale(year)[,1]) %>%
  select(year, year_scaled) %>%
  distinct(year, year_scaled) %>%
  arrange(year)

USDA_lookup <- seedf %>%
  mutate(USDA_code_id=as.integer(factor(USDA_code))) %>%
  select(USDA_code_id, USDA_code, USDA_name) %>%
  distinct() %>%
  arrange(USDA_code_id)

veg_plot_lookup <- vegf %>%
  mutate(
    plot_id=as.integer(factor(plot)),
    habitat_id=habitat=="xeric"
    ) %>%
  select(plot_id, plot, habitat_id, habitat) %>%
  distinct() %>%
  arrange(plot_id)

seed_plot_lookup <- seedf %>%
  mutate(
    plot_id=as.integer(factor(plot)),
    habitat_id=habitat=="xeric"
  ) %>%
  select(plot_id, plot, habitat_id, habitat) %>%
  distinct() %>%
  arrange(plot_id)

# Put data into list for cmdstanr model
dlist <- list(
  # Vegetation data
  Nyear=nrow(vegf),
  veg_year=vegf$year-min(vegf$year),
  Kspecies=nrow(USDA_lookup),
  Nveg=nrow(vegf),
  veg_species=as.integer(factor(vegf$USDA_code)),
  veg_count=vegf$count_scaled,
  Nveg_plots=nrow(veg_plot_lookup),
  veg_plot=as.integer(factor(vegf$plot)),
  veg_habitat=veg_plot_lookup$habitat_id,

  # Seed data
  Nseed=nrow(seedf),
  seed_change=seedf$change/(2023-1989),
  seed_species=as.integer(factor(seedf$USDA_code)),
  Nseed_plots=nrow(seed_plot_lookup),
  seed_plot=as.integer(factor(seedf$plot)),
  seed_habitat=seed_plot_lookup$habitat_id
)

# MODEL ------------------------------------------------------------------------
options(mc.cores=ifelse(parallel::detectCores()>4, 4, 2)) #set cores for parallel processing

# Compile model
mod4 <- cmdstan_model("Zenodo_archiving/scripts/stan/compare_seed_veg_change_w_habitat.stan")

# Fit model
## will get warnings as model starts sampling at extreme values but fit is fine
fit4 <- mod4$sample(
  data=dlist,
  chains=4,
  seed=5336,
  adapt_delta=.95
)

# Check model diagnostics (model doesn't always automatically warn you when there are issues)
fit4$cmdstan_diagnose()

# FIGURES ----------------------------------------------------------------------
theme_set(ggthemes::theme_tufte() + theme(panel.border=element_rect(fill=NA))) #set theme for plotting

## find rank order for change in seed abundance to sort plot later
sp_order <- lapply(USDA_lookup$USDA_code_id, function(x) {
  tidy_draws(fit4) %>%
    mutate(
      USDA_code_id=x,
      # Expected seed change WITHOUT habitat effects (baseline)
      epred_base=get(paste0("a_seed[", x, "]")) + b_seed * get(paste0("b_veg[", x, "]")),

      # Expected seed change for XERIC habitat (habitat=0)
      epred_mesic=get(paste0("a_seed[", x, "]")) +
        b_seed * get(paste0("b_veg[", x, "]")) +
        get(paste0("beta_seed_habitat[", x, "]")) * 0 +
        get(paste0("beta_seed_habitat_bveg[", x, "]")) * 0 * get(paste0("b_veg[", x, "]")),

      # Expected seed change for MESIC habitat (habitat=1)
      epred_xeric=get(paste0("a_seed[", x, "]")) +
        (b_seed + get(paste0("beta_seed_habitat_bveg[", x, "]"))) * get(paste0("b_veg[", x, "]")) +
        get(paste0("beta_seed_habitat[", x, "]")) * 1
    ) %>%
    select(USDA_code_id, epred_base, epred_xeric, epred_mesic)
}) %>%
  do.call(rbind, .) %>%
  mutate(
    USDA_code=USDA_lookup$USDA_code[USDA_code_id],
    USDA_name=USDA_lookup$USDA_name[USDA_code_id]) %>%
  group_by(USDA_code, USDA_name) %>%
  summarize(mean=mean(epred_base), .groups="drop") %>%
  arrange(mean) %>%
  select(USDA_code, USDA_name)

## get species names for plotting
sp_code_name <- read.csv("data/seedbank_composition.ac_hh.data.csv") %>%
  select(USDA_code, USDA_name) %>% #select sp codes and names
  distinct() %>% #remove dupes
  filter(USDA_code%in%sp_order) %>%  #remove taxa that weren't in analysis
  arrange(factor(USDA_code, levels=sp_order))

# Posterior distributions for change in seed abundance for each species
figS3A <- lapply(USDA_lookup$USDA_code_id, function(x) {
  tidy_draws(fit4) %>%
    mutate(
      USDA_code_id=x,
      # Expected seed change for MESIC habitat (habitat=0 -> b_veg col 1)
      mesic=get(paste0("a_seed[", x, "]")) +
        b_seed_mesic * get(paste0("b_veg[", x, ",1]")) +
        get(paste0("beta_seed_habitat[", x, "]")) * 0 +
        get(paste0("beta_seed_habitat_bveg[", x, "]")) * 0 * get(paste0("b_veg[", x, ",1]")),
      # Expected seed change for XERIC habitat (habitat=1 -> b_veg col 2)
      xeric=get(paste0("a_seed[", x, "]")) +
        (b_seed_mesic + get(paste0("beta_seed_habitat_bveg[", x, "]"))) * get(paste0("b_veg[", x, ",2]")) +
        get(paste0("beta_seed_habitat[", x, "]")) * 1
    ) %>%
    select(USDA_code_id, mesic, xeric)
}) %>%
  do.call(rbind, .) %>%
  pivot_longer(cols=c(mesic, xeric), names_to="habitat", values_to="epred") %>%
  mutate(USDA_name=factor(USDA_lookup$USDA_name[USDA_code_id], levels=sp_order$USDA_name)) %>%
  ggplot(aes(x=epred, y=fct_rev(USDA_name))) +
  geom_hline(yintercept=sp_order$USDA_name, linetype=2, color="grey25", linewidth=0.1) +
  stat_slab(fill="black") +
  geom_vline(xintercept=0, linetype=2, color="red") +
  scale_x_continuous(name="Seed Change") +
  scale_y_discrete(name="Species") +
  facet_grid(~habitat) +
  theme(axis.text.y=element_text(face="italic"))

# Posterior distributions for change in veg cover for each species
figS3B <- fit4$draws("b_veg", format="df") %>%
  select(starts_with("b_veg")) %>%
  pivot_longer(cols=everything(), names_to="parameter", values_to="b_veg") %>%
  mutate(
    USDA_code_id=as.integer(str_extract(parameter, "(?<=\\[)\\d+")),
    USDA_name=factor(USDA_lookup$USDA_name[USDA_code_id], levels=sp_order$USDA_name),
    habitat_id=as.integer(str_extract(parameter, "(?<=,)\\d+")),
    habitat=if_else(habitat_id == 1, "mesic", "xeric")
  ) %>%
  ggplot(aes(x=b_veg, y=fct_rev(USDA_name))) +
  geom_hline(yintercept=sp_order$USDA_name, linetype=2, color="grey25", linewidth=0.1) +
  stat_slab(fill="black", normalize="groups") +
  geom_vline(xintercept=0, linetype=2, color="red") +
  scale_x_continuous(name="Veg Change") +
  facet_grid(~habitat) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

# Posterior distribution for relationship between seed and veg change
figS3C <- tidy_draws(fit4) %>%
  # compute mean xeric slope across species by averaging beta_seed_habitat_bveg[1..K]
  mutate(
    b_seed_xeric=b_seed_mesic + rowMeans(
      across(starts_with("beta_seed_habitat_bveg"))
    )
  ) %>%
  select(b_seed_mesic, b_seed_xeric) %>%
  pivot_longer(
    cols=everything(),
    names_to="habitat",
    values_to="b_seed"
  ) %>%
  mutate(habitat=recode(habitat,
                          b_seed_mesic ="mesic",
                          b_seed_xeric ="xeric"
  )) %>%
  ggplot(aes(x=b_seed)) +
  stat_slab(alpha=0.6, normalize="groups") +
  geom_vline(xintercept=0, linetype=2, color="red") +
  scale_x_continuous(name="Seed Change ~ Veg Change") +
  facet_grid(~habitat) +
  theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
    )

figS3 <- (figS3A | figS3B) / figS3C +
  plot_annotation(tag_levels='A')

# Write Figure 3
ggsave("Zenodo_archiving/figures/figS3.pdf", figS3, width=7.5, height=5, units="in", dpi=600)
