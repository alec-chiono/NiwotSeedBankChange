# Visualizations of posteriors for Bayesian model evaluating relationship between long-term above-ground vegetation changes and seed bank changes
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ----
library(tidybayes)
library(ggdist)
library(patchwork)

# MODEL ----
source("scripts/04_LongTermSeedVegChange_02_fit.R")

# POSTERIOR PREDICTIONS ----
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

# VIZ ----
## Set plotting theme
theme_set(
  ggthemes::theme_tufte() +
    theme(
      panel.border = element_rect(fill = NA),
      axis.title.x = element_text(size = 8)
    )
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
fig4A <- seed_draws %>%
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
fig4B <- veg_draws %>%
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
fig4C <- rel_draws %>%
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

fig4 <- ((fig4A + fig4B + plot_layout(axes = "collect")) / fig4C) +
  plot_layout(guides = "collect")

ggsave(
  "figures/fig4.pdf",
  fig4,
  width = 7.5,
  height = 5,
  units = "in",
  dpi = 600
)
