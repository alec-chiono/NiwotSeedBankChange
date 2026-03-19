# Compare change in seed bank community composition over time using redundancy analysis (RDA)
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ---------------------------------------------------------------------
librarian::shelf(tidyverse, vegan, permute, betapart, indicspecies, ggrepel)

# DATA -------------------------------------------------------------------------
## Download
source("scripts/source/download_data.R")

## Wrangle
sc_df <- data_list$seedbank_composition.ac_hh.data %>%
  filter(substr(USDA_code,1, 1)!=2 & USDA_code!="CAREX" & USDA_code!="POA") %>%  #remove records not identified to species
  select(year:USDA_code, count) %>% #select relevant columns
  pivot_wider(names_from=USDA_code, values_from=count) %>%  #pivot into proper format for analyses
  filter(rowSums(across(-(year:depth))) > 0) #remove samples that had no species present

# RDA ---------------------------------------------------------------------
## Pull out abundance data
sc_matrix <- sc_df %>%
  select(-(year:depth)) %>%
  as.matrix()

## Pull out meta data
sc_meta <- sc_df %>%
  select(year:depth) %>%
  mutate(across(everything(), as.factor)) #pull out just metadata and make year a factor

## DCA to determine if we should do RDA or CCA
dca <- decorana(sc_matrix)
print(dca) #DCA1 axis length is 1.003, so will do RDA

## Fit RDA
sc_rda <- rda(decostand(sc_matrix, method="hellinger") ~ #Hellinger transformation (https://r.qcbs.ca/workshop09/book-en/transformations.html)
                year*habitat*depth #variables of interest
                + Condition(plot), #spatial structure
              data=sc_meta)

## Evaluate model significance
anova.cca(sc_rda, step=1000) #check overall significance
anova.cca(sc_rda, step=1000, by="axis") #check axes significance
anova.cca(sc_rda, step=1000, by="margin", scope=c("year", "habitat", "depth", "year:habitat", "year:depth", "habitat:depth", "year:habitat:depth")) #check which terms significant

## Extract sample scores
sc_scores <- sc_meta %>%
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

# PERMANOVA --------------------------------------------------------------------
## Bray-Curtis distance matrix
bc_dist <- vegdist(sc_matrix, method="bray")

## Fit PERMANOVA
adonis2(bc_dist ~
          year*habitat*depth,
        data=sc_meta,
        permutations=how(nperm=999, blocks=sc_meta$plot),
        by="terms")

# BETADISPER -------------------------------------------------------------------
## Test for homogeneity of multivariate dispersions (i.e., beta diversity)

# Year
bd_year <- betadisper(bc_dist, group=sc_meta$year)
permutest(bd_year, permutations=999)

# Habitat
bd_habitat <- betadisper(bc_dist, group=sc_meta$habitat)
permutest(bd_habitat, permutations=999)

# Depth
bd_depth <- betadisper(bc_dist, group=sc_meta$depth)
permutest(bd_depth, permutations=999)

# BETAPART DECOMPOSITION -------------------------------------------------------
# Requires raw abundance matrix (not hellinger transformed)
beta_abund <- beta.pair.abund(sc_matrix, index.family="bray")

# Three distance matrices:
# beta_abund$beta.bray    =total Bray-Curtis dissimilarity
# beta_abund$beta.bray.bal=balanced variation (turnover component)
# beta_abund$beta.bray.gra=abundance gradient component

# Total
adonis2(beta_abund$beta.bray ~ year*habitat*depth,
        data=sc_meta, permutations=how(nperm=999, blocks=sc_meta$plot), by="terms")

# Turnover component
adonis2(beta_abund$beta.bray.bal ~ year*habitat*depth,
        data=sc_meta, permutations=how(nperm=999, blocks=sc_meta$plot), by="terms")

# Abundance gradient component
adonis2(beta_abund$beta.bray.gra ~ year*habitat*depth,
        data=sc_meta, permutations=how(nperm=999, blocks=sc_meta$plot), by="terms")



# INDICATOR SPECIES ANALYSIS ---------------------------------------------------
# Create a combined grouping variable for year x habitat combinations
sc_meta <- sc_meta %>%
  mutate(year_habitat = interaction(year, habitat, sep = "_"))

# Run indicator species analysis
# Using the abundance matrix (not presence/absence)
indval_result <- multipatt(sc_matrix,
                           cluster = sc_meta$year_habitat,
                           func = "IndVal.g",        # corrected IndVal
                           control = how(nperm=999, blocks=sc_meta$plot))

summary(indval_result, indvalcomp = TRUE, alpha = 0.05)
