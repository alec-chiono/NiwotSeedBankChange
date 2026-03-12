# Compare co-located above-ground vegetation and seed bank community composition using distance-based redundancy analysis (dbRDA)
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ---------------------------------------------------------------------
librarian::shelf(tidyverse, vegan)

# DATA -------------------------------------------------------------------------
## Download
source("Zenodo_archiving/scripts/source/download_data.R")

## Wrangle data
### Seed bank data
seed_df <- data_list$seedbank_composition.ac_hh.data %>%
  filter(
    year==2023, #use only recent data
    substr(USDA_code,1, 1)!=2 & USDA_code!="CAREX" & USDA_code!="POA" #remove records not identified to species
   ) %>%
  group_by(habitat, plot, replicate, USDA_code) %>%
  summarize(count=sum(count), .groups="drop") %>% #sum counts across depths
  mutate(
    community="seed", #denote this as seed community data
    present=ifelse(count>0, 1, 0) #convert count data into presence-absence data
    ) %>%
  select(community, habitat:USDA_code, present) #select relevant columns

### Co-located veg data
veg_df <- data_list$veg_composition.ac_hh.data %>%
  filter(substr(USDA_code,1, 1)!=2 & USDA_code!="POA" & USDA_code!="CAREX") %>% #remove records not identified to species
  mutate(community="veg") %>% #denote this as vegetation community data
  select(community, habitat:plot, USDA_code) %>% #select relevant columns
  distinct() %>% mutate(present=1) #convert into presence data

### Combine seed and veg data
sv_df <- bind_rows(seed_df, veg_df) %>%
  complete(nesting(community, habitat, plot, replicate), USDA_code, fill=list(present=0)) %>% #turn into complete presence-absence data set
  pivot_wider(names_from=USDA_code, values_from=present) #pivot into wide format for use in RDA

# dbRDA ------------------------------------------------------------------------
## Pull out community data and convert into Jaccard distance matrix
sv_dist <- as.matrix(select(sv_df, -(community:replicate))) %>%
  vegdist(method="jaccard")

## Fit dbRDA
sv_rda <- dbrda(sv_dist ~
                  community*habitat #variables of interest
                + Condition(plot), #spatial structure
                data=sv_df,
                distance="jaccard")

# Evaluate model significance
anova.cca(sv_rda, step=1000) #check overall significance
anova.cca(sv_rda, step=1000, by="axis") #check axes significance
anova.cca(sv_rda, step=1000, by="margin", scope=c("community", "habitat", "community:habitat")) #check which terms are significant

## Extract scores
sv_scores <- sv_df %>%
  select(community:replicate) %>%
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
  scale_color_manual(values=c("brown4", "green4"))

# Write Figure 2
ggsave("Zenodo_archiving/figures/fig2.pdf", fig2, width=7.5, height=5, units="in", dpi=600)
