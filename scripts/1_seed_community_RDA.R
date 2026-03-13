# Compare change in seed bank community composition over time using redundancy analysis (RDA)
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ---------------------------------------------------------------------
librarian::shelf(tidyverse, vegan, ggrepel)

# DATA -------------------------------------------------------------------------
## Download
source("scripts/source/download_data.R")

## Wrangle
sc_df <- data_list$seedbank_composition.ac_hh.data %>%
  filter(substr(USDA_code,1, 1)!=2 & USDA_code!="CAREX" & USDA_code!="POA") %>%  #remove records not identified to species
  select(year:USDA_code, count) %>% #select relevant columns
  pivot_wider(names_from=USDA_code, values_from=count) #pivot into format for RDA

# RDA --------------------------------------------------------------------------
## Pull out community data and do Hellinger transformation
sc_matrix <- as.matrix(select(sc_df, -(year:depth))) %>% #pull out just community data for use in rda()
  decostand(method = "hellinger") #and do Hellinger transformation (https://r.qcbs.ca/workshop09/book-en/transformations.html)

## Fit RDA
sc_rda <- rda(sc_matrix ~
                year*habitat*depth #variables of interest
                + Condition(plot), #spatial structure
              data=sc_df)

## Evaluate model significance
anova.cca(sc_rda, step=1000) #check overall significance
anova.cca(sc_rda, step=1000, by="axis") #check axes significance
anova.cca(sc_rda, step=1000, by="margin", scope=c("year", "habitat", "depth", "year:habitat", "year:depth", "habitat:depth", "year:habitat:depth")) #check which terms significant

## Extract sample scores
sc_scores <- sc_df %>%
  select(year:depth) %>%
  mutate(year=factor(year)) %>%
  cbind(scores(sc_rda)$sites)

## Extract species scores
spp_scores <- scores(sc_rda)$species %>%
  as.data.frame() %>%
  rownames_to_column("USDA_code")

# FIGURES ----------------------------------------------------------------------
## Set plotting theme
theme_set(ggthemes::theme_tufte() + theme(panel.border=element_rect(fill=NA)))

## Fig. 1A: Just sample scores
fig1A <- sc_scores %>%
  mutate(habitat_depth=paste(habitat, depth, sep="_")) %>%
  ggplot(aes(x=RDA1, y=RDA2)) +
  geom_vline(xintercept=0, linetype="dashed", color="gray70") +
  geom_hline(yintercept=0, linetype="dashed", color="gray70") +
  geom_point(aes(color=year, shape=habitat), size=5) +
  scale_shape_manual(values=c(16, 17)) +
  scale_color_manual(values=c("skyblue3", "red4"))

## Fig. S1: Sample scores and species scores
figS1 <- sc_scores %>%
  mutate(habitat_depth=paste(habitat, depth, sep="_")) %>%
  ggplot(aes(x=RDA1, y=RDA2)) +
  geom_vline(xintercept=0, linetype="dashed", color="gray70") +
  geom_hline(yintercept=0, linetype="dashed", color="gray70") +
  geom_point(aes(color=year, shape=habitat), size=5) +
  scale_shape_manual(values=c(16, 17)) +
  scale_color_manual(values=c("skyblue3", "red4")) +
  geom_label_repel(data=spp_scores, aes(label=USDA_code), size=3, max.overlaps=100)

### Write Fig. S1
ggsave("figures/figS1.pdf", figS1, width=7.5, height=7.5, units="in", dpi=600)

# Distance-based RDA -----------------------------------------------------------
## to evaluate if compositional changes are due to species turnover

## Wrangle data for dbRDA
sc_dist_df <- sc_df %>%
  mutate(row_sum = rowSums(select(., -(year:depth)), na.rm = TRUE)) %>%
  filter(row_sum > 0) %>% #remove samples that did not have any species present (since these would have distance of 0 to all other samples and cause issues for dbRDA)
  select(-row_sum) %>%
  mutate(across(-(year:depth), ~ ifelse(. > 0, 1, 0))) #convert to presence absence

sc_dist <- sc_dist_df %>%
  select(-(year:depth)) %>%
  as.matrix() %>%
  vegdist(method="jaccard") #convert into jaccard distance matrix

## Fit dbRDA
sc_dbrda <- dbrda(sc_dist ~
                    year*habitat*depth #variables of interest
                  + Condition(plot), #spatial structure
                  data=sc_dist_df,
                  distance="jaccard")

## Evaluate model significance
anova.cca(sc_dbrda, step=1000) #check overall significance
anova.cca(sc_dbrda, step=1000, by="axis") #check axes significance
anova.cca(sc_dbrda, step=1000, by="margin", scope=c("year", "habitat", "depth", "year:habitat", "year:depth", "habitat:depth", "year:habitat:depth")) #check which terms significant

## Extract sample scores
sc_db_scores <- sc_dist_df %>%
  select(year:depth) %>%
  cbind(scores(sc_dbrda)$sites) %>%
  mutate(year=factor(year))

## Fig. S2: Species composition based on presence-absence alone
figS2 <- sc_db_scores %>%
  mutate(habitat_depth=paste(habitat, depth, sep="_")) %>%
  ggplot(aes(x=dbRDA1, y=dbRDA2 , group=paste0(habitat,plot,depth))) +
  geom_vline(xintercept=0, linetype="dashed", color="gray70") +
  geom_hline(yintercept=0, linetype="dashed", color="gray70") +
  geom_point(aes(color=year, shape=habitat), size=5) +
  scale_shape_manual(values=c(16, 17)) +
  scale_color_manual(values=c("skyblue3", "red4"))

### Write Fig. S2
ggsave("figures/figS2.pdf", figS2, width=7.5, height=7.5, units="in", dpi=600)
