// Heirarchical negative binomial model to estimate species counts in seed bank
// and then derive biodiversity indices: richness, evenness, Hill N1
// Alec Chiono; alec.chiono@colorado.edu

data {
  int<lower=1> N_obs;                    // total number of observations
  int<lower=1> N_species;                // number of species
  int<lower=1> N_plots;                  // number of sampling plots
  int<lower=1> N_years;                  // number of years in the study
  int<lower=1> N_habitats;               // number of habitat types
  array[N_obs] int<lower=0> count;       // observed count (response variable)
  array[N_obs] int<lower=1> species;     // species ID for each observation
  array[N_obs] int<lower=1> plot;        // plot ID for each observation
  array[N_obs] int<lower=1> year_idx;    // year index for each observation
  array[N_obs] int<lower=1> habitat_idx; // habitat index for each observation
}

parameters {
  real mu_global;  // global intercept
  real year_effect;        // fixed effect for year (2023 vs. 1989 contrast)
  real habitat_effect;     // fixed effect for habitat (xeric vs. mesic contrast)
  real year_hab_int;       // fixed year × habitat interaction
  vector[N_plots] plot_effect_raw;             // standardized plot effects
  real<lower=0> sigma_plot;                    // SD of plot effects
  matrix[N_species, N_habitats] species_habitat_raw;        // baseline habitat:species interaction
  matrix[N_species, N_habitats] species_habitat_change_raw; // change in habitat:species interaction across years
  real<lower=0> sigma_species_habitat;         // SD of baseline habitat:species interaction
  real<lower=0> sigma_species_habitat_change;  // SD of change in habitat:species interaction across years
  real<lower=0> phi;
}

transformed parameters {
  // Linear predictor (log-scale abundance/intensity)
  vector[N_obs] log_lambda;
  // Apply SD to mean effects to generate actual effects
  vector[N_plots] plot_effect = sigma_plot * plot_effect_raw;
  matrix[N_species, N_habitats] species_habitat = sigma_species_habitat * species_habitat_raw;
  matrix[N_species, N_habitats] species_habitat_change = sigma_species_habitat_change * species_habitat_change_raw;


  // Construct the expected log-abundance for each observation
  for (i in 1:N_obs) {
    // Temporal change: apply only to the most recent year
    real temporal_change = (year_idx[i] == N_years) ? species_habitat_change[species[i], habitat_idx[i]] : 0;

    log_lambda[i] =
      mu_global +
      year_effect * (year_idx[i] == N_years ? 1 : 0) +
      habitat_effect * (habitat_idx[i] == 2 ? 1 : 0) +
      year_hab_int * (year_idx[i] == N_years ? 1 : 0) * (habitat_idx[i] == 2 ? 1 : 0) +
      plot_effect[plot[i]] +
      species_habitat[species[i], habitat_idx[i]] +
      temporal_change;
  }
}

model {
  // Priors (weakly informative)
  // Global intercept
  mu_global ~ student_t(3, 0, 2.5);
  // Fixed effects
  year_effect ~ student_t(3, 0, 2.5);
  habitat_effect ~ student_t(3, 0, 2.5);
  year_hab_int ~ student_t(3, 0, 1);
  // SDs slightly tighter
  sigma_plot ~ normal(0, 0.5);
  sigma_species_habitat ~ normal(0, 0.5);
  sigma_species_habitat_change ~ normal(0, 0.25); // interactions even tighter
  // All "raw" terms standard normal
  plot_effect_raw ~ std_normal();
  to_vector(species_habitat_raw) ~ std_normal();
  to_vector(species_habitat_change_raw) ~ std_normal();
  // Overdispersion
  phi ~ exponential(1);


  // Likelihood
  count ~ neg_binomial_2_log(log_lambda, phi);
}

generated quantities {
  vector[N_obs] log_lambda_gq = log_lambda;
  vector[N_obs] log_lik;
  for (n in 1:N_obs) log_lik[n] = neg_binomial_2_log_lpmf(count[n] | log_lambda[n], phi);

  array[N_years, N_habitats] real richness;
  array[N_years, N_habitats] real evenness;
  array[N_years, N_habitats] real hill_N1;
  array[N_years, N_habitats, N_species] int predicted_seed_counts;

  for (y in 1:N_years) {
    for (h in 1:N_habitats) {

      // Log-scale expected abundance per species
      vector[N_species] log_expected =
        mu_global +
        year_effect * (y == N_years ? 1 : 0) +
        habitat_effect * (h == 2 ? 1 : 0) +
        year_hab_int * (y == N_years ? 1 : 0) * (h == 2 ? 1 : 0) +
        species_habitat[, h] +
        (y == N_years ? species_habitat_change[, h] : rep_vector(0, N_species));

      // ── Step 1: Draw predicted counts FIRST ──────────────────────────
      for (s in 1:N_species)
        predicted_seed_counts[y, h, s] = neg_binomial_2_rng(exp(log_expected[s]), phi);

      // ── Step 2: Diversity metrics conditioned on predicted counts ≥ 1 ─
      // Sum predicted counts across present species for proportion calculation
      // In the y/h loop, after computing log_expected:
      vector[N_species] lambda_s = exp(log_expected);
      vector[N_species] p_present;  // P(count >= 1) per species

      for (s in 1:N_species)
        p_present[s] = 1.0 - exp(neg_binomial_2_lpmf(0 | lambda_s[s], phi));

      // Expected richness: sum of presence probabilities (smooth, continuous)
      real richness_real = sum(p_present);
      richness[y, h] = richness_real;

      // Soft-weighted abundance: down-weight species unlikely to be present
      vector[N_species] weighted_abund = p_present .* lambda_s;
      real total_weighted = sum(weighted_abund);

      if (total_weighted > 0) {
        real thresh = 0.1;  // species present in >10% of hypothetical samples

        real total_thresh = 0.0;
        int S_int = 0;
        real H_present = 0.0;

        for (s in 1:N_species)
          if (p_present[s] > thresh)
            total_thresh += weighted_abund[s];

        for (s in 1:N_species) {
          if (p_present[s] > thresh) {
            S_int += 1;
            if (total_thresh > 0) {
              real w = weighted_abund[s] / total_thresh;
              if (w > 1e-12)
                H_present += -w * log(w);
            }
          }
        }

        // Richness: soft expected value for all species (keep as before)
        richness[y, h] = richness_real;

        // H_max and evenness now use the same S as H_present
        real H_max = (S_int > 1) ? log(S_int) : 0.0;
        evenness[y, h] = (H_max > 0) ? H_present / H_max : 0.0;

        // Hill N1 can stay soft-weighted (no H_max involved, so no overflow risk)
        real H_hill = 0.0;
        for (s in 1:N_species) {
          real w = weighted_abund[s] / total_weighted;
          if (w > 1e-12)
            H_hill += -w * log(w);
        }
        hill_N1[y, h] = exp(H_hill);

      } else {
        richness[y, h] = 0.0;
        evenness[y, h] = 0.0;
        hill_N1[y, h] = 0.0;
      }
    }
  }
}
