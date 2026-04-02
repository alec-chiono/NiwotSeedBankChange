# Compare co-located above-ground vegetation and seed bank community composition using distance-based redundancy analysis (dbRDA)
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ---------------------------------------------------------------------
librarian::shelf(tidyverse, vegan, indicspecies)

# DATA -------------------------------------------------------------------------
## Download
source("scripts/source/download_data.R")

## Wrangle data
### Seed bank data
seed_df <- data_list$seedbank_composition.ac_hh.data %>%
  filter(
    year==2023, #use only recent data
    substr(USDA_code,1, 1)!=2 & USDA_code!="CAREX" & USDA_code!="POA" #remove records not identified to species
   ) %>%
  group_by(habitat, plot, replicate, USDA_name) %>%
  summarize(count=sum(count), .groups="drop") %>% #sum counts across depths
  mutate(
    community="seed", #denote this as seed community data
    present=ifelse(count>0, 1, 0) #convert count data into presence-absence data
    ) %>%
  select(community, habitat:USDA_name, present) #select relevant columns

### Co-located veg data
veg_df <- data_list$veg_composition.ac_hh.data %>%
  filter(substr(USDA_code,1, 1)!=2 & USDA_code!="POA" & USDA_code!="CAREX") %>% #remove records not identified to species
  mutate(community="vegetation") %>% #denote this as vegetation community data
  select(community, habitat:plot, USDA_name) %>% #select relevant columns
  distinct() %>% mutate(present=1) #convert into presence data

### Combine seed and veg data
sv_df <- bind_rows(seed_df, veg_df) %>%
  complete(nesting(community, habitat, plot, replicate), USDA_name, fill=list(present=0)) %>% #turn into complete presence-absence data set
  pivot_wider(names_from=USDA_name, values_from=present) %>% #pivot into wide format for use in RDA
  mutate(across(community:replicate, as.factor)) #make metadata columns factors

# dbRDA ------------------------------------------------------------------------
## Pull out presence/absence data and convert into Jaccard distance matrix
sv_dist <- sv_df %>%
  select(-(community:replicate)) %>%
  as.matrix() %>%
  vegdist(method="jaccard")

## Pull out meta data
sv_meta <- sv_df %>%
  select((community:replicate)) %>%
  mutate(across(everything(), as.factor))

## Fit dbRDA
sv_rda <- dbrda(sv_dist ~
                  community*habitat #variables of interest
                + Condition(plot), #spatial structure
                data=sv_meta,
                distance="jaccard")

# Evaluate model significance
anova.cca(sv_rda, step=1000) #check overall significance
anova.cca(sv_rda, step=1000, by="axis") #check axes significance
anova.cca(sv_rda, step=1000, by="margin", scope=c("community", "habitat", "community:habitat")) #check which terms are significant

## Extract scores
sv_scores <- sv_meta %>%
  cbind(scores(sv_rda)$sites)

# FIGURES ----------------------------------------------------------------------
## Set plotting theme
theme_set(ggthemes::theme_tufte() + theme(panel.border=element_rect(fill=NA)))

## Fig. 2: dbRDA comparing seed bank and veg community from same locations
fig2 <- sv_scores %>%
  ggplot(aes(x=dbRDA1, y=dbRDA2)) +
  geom_vline(xintercept=0, linetype="dashed", color="gray70") +
  geom_hline(yintercept=0, linetype="dashed", color="gray70") +
  geom_point(aes(color=community, shape=habitat), size=5) +
  scale_shape_manual(values=c(16, 17)) +
  scale_color_manual(values=c("brown4", "green4")) +
  theme(legend.position="top")

# Write Figure 2
ggsave("figures/fig2.pdf", fig2, width=7.5, height=5, units="in", dpi=600, bg="white")

# PERMANOVA --------------------------------------------------------------------
## Aggregate soil samples for each plot so blocking design does not invalidate PERMANOVA
sv_df_agg <- sv_df %>%
  group_by(community, habitat, plot) %>%
  summarize(across(where(is.numeric), max), .groups = "drop")
# max() preserves presence if species present in either replicate

## Rebuild distance matrix and meta df
sv_dist_agg <- sv_df_agg %>%
  select(-(community:plot)) %>%
  as.matrix() %>%
  vegdist(method="jaccard")

sv_meta_agg <- sv_df_agg %>%
  select((community:plot)) %>%
  mutate(across(everything(), as.factor))


## PERMANOVA on aggregated data
adonis2(sv_dist_agg ~ community * habitat,
        data = sv_meta_agg,
        permutations = how(blocks = sv_meta_agg$plot, nperm = 999),
        by = "terms")

## Dispersion tests
bd_community <- betadisper(sv_dist_agg, group = sv_meta_agg$community)
permutest(bd_community, permutations = 999)

bd_habitat <- betadisper(sv_dist_agg, group = sv_meta_agg$habitat)
permutest(bd_habitat, permutations = 999)

# INDICATOR SPECIES ANALYSIS ---------------------------------------------------
# Use aggregated data and community type as grouping variable
sv_indval_combined <- multipatt(
  select(sv_df_agg, -(community:plot)),  # adjust exclusion range as needed
  cluster = sv_df_agg$community,
  func = "IndVal.g",
  control = how(nperm = 999)
)

summary(sv_indval_combined, indvalcomp = TRUE, alpha = 0.05)

sv_indval_formatted <- sv_indval_combined$sign %>%
  rownames_to_column("species_code") %>%
  filter(p.value <= 0.05) %>%
  mutate(
    group = case_when(
      s.seed       == 1 & s.vegetation == 0 ~ "Seed bank",
      s.vegetation == 1 & s.seed       == 0 ~ "Vegetation",
      TRUE ~ "Both"
    ),
    stat    = round(stat, 3),
    p.value = round(p.value, 3)
  ) %>%
  select(species_code, group, stat, p.value) %>%
  arrange(group, desc(stat)) %>%
  mutate(p.value = ifelse(p.value <= 0.001, "< 0.001", as.character(p.value)))

ft_sv <- flextable(sv_indval_formatted) %>%
  set_header_labels(
    values = setNames(
      c("Species", "Community", "IndVal.g", "p"),
      names(sv_indval_formatted)
    )
  ) %>%
  italic(j = "species_code") %>%
  bold(part = "header") %>%
  hline(i = which(diff(as.numeric(factor(sv_indval_formatted$group))) != 0)) %>%
  align(j = c("stat", "p.value"), align = "center", part = "all") %>%
  align(j = c("species_code", "group"), align = "left", part = "all") %>%
  set_caption("Table S4. Indicator species associated with seed bank and vegetation communities (IndVal.g; Cáceres & Legendre 2009). Only species with p \u2264 0.05 shown.") %>%
  theme_booktabs() %>%
  autofit()

save_as_docx(ft_sv, path = "tables/tableS4.docx")
