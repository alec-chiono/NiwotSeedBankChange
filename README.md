# NiwotSeedBankChange

Code and data for analyses of long-term seed bank and vegetation changes in the Niwot Ridge alpine tundra ecosystem.

**Author:** Alec Chiono  
**Affiliation:** Department of Ecology and Evolutionary Biology, University of Colorado Boulder  
**Contact:** alec.chiono@colorado.edu  
**License:** The Unlicense (public domain)

---

## Associated publication

> Chiono, A., Humphries, H. C., Willis, G., & Emery, N. C. (in prep). Seed bank reorganization in an alpine dry meadow plant community.

---

## Repository structure

```
NiwotSeedBankChange/
├── data/                                    # Input data files and metadata
│   └── saddptqd_xericmesic_categorization.csv  # Habitat classification (mesic/xeric) for Saddle Grid plots
│
├── scripts/                                 # R analysis scripts (organized by analysis module)
│   ├── 00_Function_data_download.R         # Helper: download data from EDI repository
│   ├── 00_Function_plot_prc_hist.R         # Helper: posterior retrodictive check plots & parameter summaries
│   │
│   ├── 01_Figure1.R                        # Descriptive: bar plot of seed counts by species & year (Figure 1)
│   │
│   ├─── MODULE 2: Seed Bank Community Change (1989 vs 2023)
│   ├── 02_SeedCommunityChange_01_prep.R        # Prep seed bank data for multivariate analyses
│   ├── 02_SeedCommunityChange_02_PERMANOVA.R   # PERMANOVA: test for compositional differences
│   ├── 02_SeedCommunityChange_03_BETAPART.R    # Beta diversity: decompose turnover vs abundance shift
│   ├── 02_SeedCommunityChange_04_RDA.R         # RDA: ordinate community composition by year × habitat (Figure 2A)
│   │
│   ├─── MODULE 2 (continued): Bayesian Richness/Diversity Modeling
│   ├── 02_SeedCommunityChange_05_STAN_01_prep.R    # Prep data for Stan model
│   ├── 02_SeedCommunityChange_05_STAN_02_fit.R     # Fit negative binomial hierarchical model
│   ├── 02_SeedCommunityChange_05_STAN_03_checks.R  # Model diagnostics & convergence checks
│   ├── 02_SeedCommunityChange_05_STAN_04_viz.R     # Visualize richness, evenness, diversity (Figure 2B, 2C)
│   │
│   ├─── MODULE 3: Seed Bank vs Above-Ground Vegetation Composition
│   ├── 03_CompareColocatedSeedVeg_01_prep.R    # Align seed & veg data at colocated 2023 plots
│   ├── 03_CompareColocatedSeedVeg_02_PERMANOVA.R   # PERMANOVA: seed community ≠ veg community?
│   ├── 03_CompareColocatedSeedVeg_03_dbRDA.R      # dbRDA: ordinate seed-veg differences (Figure 3)
│   │
│   ├─── MODULE 4: Long-Term Seed–Vegetation Correlation (1989–2023)
│   ├── 04_LongTermSeedVegChange_01_prep.R      # Prep species-level change data: seed & veg trends
│   ├── 04_LongTermSeedVegChange_02_fit.R       # Fit Stan model: seed change ~ veg change × habitat
│   ├── 04_LongTermSeedVegChange_03_checks.R    # Model diagnostics & convergence checks
│   ├── 04_LongTermSeedVegChange_04_viz.R       # Main figures: species slopes & seed-veg correlation (Figure 4)
│   └── 04_LongTermSeedVegChange_05_wo_habitat.R # Sensitivity: model without habitat main effect
│
├── stan/                                    # Stan model definitions
│   ├── community_indices.stan               # Neg binomial hierarchical model: richness, evenness, diversity diversity
│   ├── longterm_change.stan                 # Linked linear regressions: seed change ~ veg change × habitat
│   └── longterm_change_wo_habitat.stan      # Sensitivity: without habitat main effect
│
├── figures/                                 # Publication-ready figures (auto-generated)
│   │
│   ├── Figure1.pdf                         # Seed counts by species (1989 vs 2023)
│   │
│   ├── Figure2.pdf                         # Community composition and diversity change (1989 vs 2023)
│   │
│   ├── Figure3.pdf                         # Seed bank vs above-ground vegetation composition (2023)
│   │
│   ├── Figure4.pdf                         # Long-term seed–vegetation change relationship (1989–2023)
│   │
│   ├── FigureS1.pdf                        # Posterior retrodictive check (PRC) for community_indices.stan
bands).
│   │
│   ├── FigureS2.pdf                        # Model validation and diagnostics for longterm_change.stan
│   │
│   ├── FigureS3.pdf                        # Sensitivity analysis: longterm_change_wo_habitat.stan
│   │
│   └── FigureS4.pdf                        # Beta diversity partitioning (betapart analysis)
│
├── images/                                  # Icons for Figure 1
│   ├── graminoid.png
│   ├── cushion.png
│   └── forb.png
│
├── NiwotSeedBankChange.Rproj
├── LICENSE
└── README.md
```

---

## Data sources

**Data are downloaded directly from the Environmental Data Initiative (EDI) LTER data repository** using the `download_data()` function in `scripts/00_Function_data_download.R`. No manual data downloads are required; the function handles authentication and file retrieval automatically.

## Analysis workflow

### Module 1: Descriptive visualizations
- **01_Figure1.R** – Bar chart showing seed counts for each species in 1989 vs 2023, with growth form icons (graminoid, cushion, forb) below species names. 

### Module 2: Seed bank community change (1989 vs 2023)
Evaluates whether the seed bank community changed significantly over 34 years and characterizes the type of change.

1. **02_SeedCommunityChange_01_prep.R** – Converts raw seed counts to a species-by-sample matrix with metadata (year, habitat, plot, depth).

2. **02_SeedCommunityChange_02_PERMANOVA.R** – Tests for significant differences in seed bank composition between 1989 and 2023 using permutation-based multivariate analysis of variance (PERMANOVA; `adonis2` function). Tests main effects and interactions of year, habitat, and depth with permutations blocked by plot.

3. **02_SeedCommunityChange_03_BETAPART.R** – Decomposes beta diversity into turnover (species replacement) and abundance gradient (nestedness) components. PERMANOVA on each component separately to determine whether compositional change is driven by species turning over vs. shifts in relative abundances of shared species.

4. **02_SeedCommunityChange_04_RDA.R** – Redundancy analysis (RDA; Hellinger-transformed data) to ordinate and visualize community composition in 1989 vs 2023. Evaluates whether habitats diverged over time and describes the multivariate structure of change.

5. **02_SeedCommunityChange_05_STAN_*** – Bayesian hierarchical modeling of species richness, Shannon evenness, and Hill's N1 diversity (effective number of equally-abundant species):
   - **_01_prep.R** – Aggregates counts by year/habitat/plot and prepares data list for Stan model
   - **_02_fit.R** – Fits negative binomial hierarchical model using cmdstanr with plot random effects and species-habitat interactions
   - **_03_checks.R** – Diagnostic checks: Rhat, effective sample size, divergent transitions, posterior retrodictive checks
   - **_04_viz.R** – Plots posterior distributions of richness, evenness, and diversity by year × habitat (Figure 2B, C)

### Module 3: Seed bank vs vegetation at colocated 2023 plots
Compares seed bank and above-ground vegetation composition at plots sampled for both datasets in 2023.

1. **03_CompareColocatedSeedVeg_01_prep.R** – Subsets seed and vegetation data to seed bank plots, converts to presence/absence for comparable scale, aligns spatial locations.

2. **03_CompareColocatedSeedVeg_02_PERMANOVA.R** – PERMANOVA on Jaccard dissimilarity to test whether seed and vegetation communities are compositionally different.

3. **03_CompareColocatedSeedVeg_03_dbRDA.R** – Distance-based RDA (dbRDA) on Jaccard dissimilarity to ordinate and visualize seed–vegetation separation (Figure 3). Examines whether seed community composition can be predicted from vegetation composition.

### Module 4: Long-term species-level seed-vegetation relationships
Evaluates whether observed changes in individual species' seed bank abundances correlate with changes in above-ground vegetation abundance (1989–2023).

1. **04_LongTermSeedVegChange_01_prep.R** – Prepares species-level change data

2. **04_LongTermSeedVegChange_02_fit.R** – Fits Bayesian linked linear regression models using cmdstanr

3. **04_LongTermSeedVegChange_03_checks.R** – Model diagnostics and posterior predictive checks.

4. **04_LongTermSeedVegChange_04_viz.R** – Plots of posteriors 

5. **04_LongTermSeedVegChange_05_wo_habitat.R** – Sensitivity analysis: refits model without habitat main effect to test robustness of results.

---

## Stan models

### `community_indices.stan`
Hierarchical negative binomial model for seed bank community composition.

### `longterm_change.stan`
Linked linear regression models for long-term species-level seed–vegetation relationship.


### `longterm_change_wo_habitat.stan`
Sensitivity analysis version of `longterm_change.stan` that removes the habitat main effect to test whether the overall seed–vegetation relationship is robust across habitats.

---

## Dependencies

All analyses require R (version 4.3+). Key packages:

```r
# Data acquisition
install.packages("EDIutils")       # Download LTER data from EDI

# Core data wrangling and visualization
install.packages("tidyverse")      # ggplot2, dplyr, tidyr, etc.
install.packages("ggthemes")       # Publication themes
install.packages("patchwork")      # Combining plots
install.packages("ggbreak")        # Axis breaks for bar plots
install.packages("ggimage")        # Embedding images in plots

# Multivariate community analysis
install.packages("vegan")          # PERMANOVA, RDA, betadiversity
install.packages("betapart")

# Bayesian modeling
install.packages("cmdstanr")       # Stan interface (install from source)

# Model diagnostics and visualization
install.packages("bayesplot")      # Posterior checks and diagnostics
```

### R version
- Analyses conducted with R 4.5.2

### Stan compiler
- cmdstanr requires a working C++ compiler. See: https://mc-stan.org/cmdstanr/

---

## Reproducibility

To reproduce all analyses from scratch:

1. **Clone or download** the repository to your local machine.

2. **Open the R project:**
   ```r
   # In RStudio, open NiwotSeedBankChange.Rproj
   # This sets the working directory to the project root
   ```

3. **Install dependencies** (see above).

4. **Run scripts in order within each module:**

   **Module 1 (descriptive):**
   ```r
   source("scripts/01_Figure1.R")  # Generates Figure 1
   ```

   **Module 2 (community change):**
   ```r
   source("scripts/02_SeedCommunityChange_01_prep.R")
   source("scripts/02_SeedCommunityChange_02_PERMANOVA.R")
   source("scripts/02_SeedCommunityChange_03_BETAPART.R")
   source("scripts/02_SeedCommunityChange_04_RDA.R")
   source("scripts/02_SeedCommunityChange_05_STAN_01_prep.R")
   source("scripts/02_SeedCommunityChange_05_STAN_02_fit.R")  # ~2–5 min
   source("scripts/02_SeedCommunityChange_05_STAN_03_checks.R")
   source("scripts/02_SeedCommunityChange_05_STAN_04_viz.R")
   ```

   **Module 3 (colocated comparison):**
   ```r
   source("scripts/03_CompareColocatedSeedVeg_01_prep.R")
   source("scripts/03_CompareColocatedSeedVeg_02_PERMANOVA.R")
   source("scripts/03_CompareColocatedSeedVeg_03_dbRDA.R")
   ```

   **Module 4 (long-term seed-veg relationships):**
   ```r
   source("scripts/04_LongTermSeedVegChange_01_prep.R")
   source("scripts/04_LongTermSeedVegChange_02_fit.R")  # ~3–10 min
   source("scripts/04_LongTermSeedVegChange_03_checks.R")
   source("scripts/04_LongTermSeedVegChange_04_viz.R")
   source("scripts/04_LongTermSeedVegChange_05_wo_habitat.R")  # ~3–10 min
   ```

5. **Output locations:**
   - Figures: `figures/`

### Notes on reproducibility

- **Data downloads:** Scripts automatically download data from EDI on first run. Internet connection required.
- **Stan compilation:** First time a Stan model is run, C++ code is compiled (~1 min). Subsequent runs use cached compiled models.
- **Computation time:** Bayesian models take 3–10 minutes to fit on a standard laptop (4+ cores, 8+ GB RAM).
- **Random seed:** Stan sampling is stochastic; posterior samples may vary slightly between runs (negligible for inference). Specific seeds can be set in model fitting scripts.
- **Dependencies on prior scripts:** Later scripts in a module depend on data prepared by earlier scripts. Ensure you run scripts in numerical order within each module.

---

## Citation

If you use this code or data, please cite:

**Chiono, A., Humphries, H. C., Willis, G., & Emery, N. C. (in prep).** Seed bank reorganization in an alpine dry meadow plant community. 

And cite the Niwot Ridge LTER datasets:

**Chiono, A., & Humphries, H. C. (2026).** Dry meadow seedbank composition, Saddle, 1989–2023. Environmental Data Initiative. https://doi.org/10.6073/pasta/[PACKAGE_ID]

**Walker, M. D., Humphries, H. C., & Niwot Ridge LTER. (2025).** Plant species composition data for Saddle grid, 1989–ongoing. Environmental Data Initiative. https://doi.org/10.6073/pasta/[PACKAGE_ID]

---

## Contact & support

For questions about the code, methods, or data in this repository, contact Alec Chiono at alec.chiono@colorado.edu.

For issues or questions about the Niwot Ridge LTER data, contact the Niwot Ridge LTER Information Manager (https://nwt.lternet.edu/our-team)

---

