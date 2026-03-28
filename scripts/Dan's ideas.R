
librarian::shelf(tidyverse, ggrepel)

# DATA -------------------------------------------------------------------------
## Download
source("scripts/source/download_data.R")

theme_set(ggthemes::theme_tufte() + theme(panel.border=element_rect(fill=NA)))

## Plot
data_list$seedbank_composition.ac_hh.data %>%
  filter(substr(USDA_code,1, 1)!=2 & USDA_code!="CAREX" & USDA_code!="POA") %>%  #remove records not identified to species
  group_by(year, USDA_name) %>%
  summarize(sum=sum(count), .groups="drop") %>%
  pivot_wider(id_cols=USDA_name, names_from=year, values_from=sum) %>%
  ggplot(aes(x=`1989`, y=`2023`, color=USDA_name)) +
  geom_hline(yintercept=0, linetype=2, color="grey") +
  geom_vline(xintercept=0, linetype=2, color="grey") +
  geom_point() +
  geom_text_repel(aes(label=USDA_name), size=3, max.overlaps=100, box.padding=.5, force=100, max.time=20) +
  scale_x_continuous(name="Seeds in 1989", limits=c(-60, NA), breaks=c(0, 25, 50, 75)) +
  scale_y_continuous(name="Seeds in 2023", limits=c(-10, NA), breaks=c(0, 25, 50, 75, 100)) +
  theme(legend.position="none")


## Plot by habitat
data_list$seedbank_composition.ac_hh.data %>%
  filter(substr(USDA_code,1, 1)!=2 & USDA_code!="CAREX" & USDA_code!="POA") %>%  #remove records not identified to species
  group_by(year, habitat, USDA_name) %>%
  summarize(sum=sum(count), .groups="drop") %>%
  pivot_wider(id_cols=c(USDA_name, habitat), names_from=year, values_from=sum) %>%
  ggplot(aes(x=`1989`, y=`2023`, color=USDA_name)) +
  geom_hline(yintercept=0, linetype=2, color="grey") +
  geom_vline(xintercept=0, linetype=2, color="grey") +
  geom_point() +
  geom_text_repel(aes(label=USDA_name), size=3, max.overlaps=100, box.padding=.5, force=100, max.time=20) +
  scale_x_continuous(name="Seeds in 1989", limits=c(-60, NA), breaks=c(0, 25, 50, 75)) +
  scale_y_continuous(name="Seeds in 2023", limits=c(-10, NA), breaks=c(0, 25, 50, 75, 100)) +
  facet_grid(~habitat) +
  theme(legend.position="none")

## Dans Fig 4 version
left_join(seed_draws, veg_draws, by=join_by(.draw, habitat, USDA_name)) %>%
  rename(seed_change=value.x, veg_change=value.y) %>%
  mutate(Species=sub("\\s*var\\..*$", "", USDA_name)) %>%
  ggplot(aes(x=veg_change, y=seed_change, color=habitat)) +
  geom_hline(yintercept=0, linetype=2, color="grey") +
  geom_vline(xintercept=0, linetype=2, color="grey") +
  geom_point(alpha=.1) +
  scale_color_manual(values=c("mesic"="green4", "xeric"="tan4"), name="Habitat") +
  facet_wrap(~Species, nrow=3) +
  theme(strip.text = element_text(face = "italic"))
