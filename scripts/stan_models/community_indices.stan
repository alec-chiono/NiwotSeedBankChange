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
  real year_hab_int;       // fixed year x habitat interaction
  vector[N_plots] plot_effect_raw;             // standardized plot effects
  real<lower=0> sigma_plot;                    // SD of plot effects
  matrix[N_species, N_habitats] species_habitat_raw;        // baseline habitat:species interaction
  matrix[N_species, N_habitats] species_habitat_change_raw; // change in habitat:species interaction across years
  real<lower=0> sigma_species_habitat;         // SD of baseline habitat:species interaction
  real<lower=0> sigma_species_habitat_change;  // SD of change in habitat:species interaction across years
  // Species-specific overdispersion (non-centered parameterization)
  vector[N_species] phi_raw;               // standardized species-level overdispersion
  real mu_log_phi;                         // mean log overdispersion across species
  real<lower=0> sigma_log_phi;             // SD of log overdispersion across species
}

transformed parameters {
  // Linear predictor (log-scale abundance/intensity)
  vector[N_obs] log_lambda;
  // Apply SD to mean effects to generate actual effects
  vector[N_plots] plot_effect = sigma_plot * plot_effect_raw;
  matrix[N_species, N_habitats] species_habitat = sigma_species_habitat * species_habitat_raw;
  matrix[N_species, N_habitats] species_habitat_change = sigma_species_habitat_change * species_habitat_change_raw;
  // Species-specific overdispersion on natural scale
  vector<lower=0>[N_species] phi = exp(mu_log_phi + sigma_log_phi * phi_raw);

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
  // Species-specific overdispersion priors
  phi_raw ~ std_normal();
  mu_log_phi ~ normal(0, 1);       // weakly informative prior on mean log-phi
  sigma_log_phi ~ normal(0, 0.5);  // partial-pooling SD across species

  // Likelihood (species-specific phi indexed by species[i])
  for (i in 1:N_obs)
    count[i] ~ neg_binomial_2_log(log_lambda[i], phi[species[i]]);
}


generated quantities {
  // --- Pointwise log-likelihood for LOO/WAIC ---
  vector[N_obs] log_lik;
  for (n in 1:N_obs)
    log_lik[n] = neg_binomial_2_log_lpmf(count[n] | log_lambda[n], phi[species[n]]);

  // --- Posterior predictive draws (year x habitat x species) ---
  // Plot effects are marginalized out; these represent community-level predictions.
  array[N_years, N_habitats, N_species] int predicted_seed_counts;

  // --- Biodiversity indices (year x habitat) ---
  array[N_years, N_habitats] real richness;   // expected species richness
  array[N_years, N_habitats] real hill_N1;    // Hill number N1 (Shannon effective species)
  array[N_years, N_habitats] real evenness;   // Pielou's J

  for (y in 1:N_years) {
    for (h in 1:N_habitats) {

      // Expected log-abundance per species, marginalizing over plot effects.
      // Plot effects average to zero by construction (mean-zero normal prior),
      // so omitting them here gives the marginal community-level expectation.
      vector[N_species] log_expected =
        mu_global +
        year_effect    * (y == N_years ? 1 : 0) +
        habitat_effect * (h == 2       ? 1 : 0) +
        year_hab_int   * (y == N_years ? 1 : 0) * (h == 2 ? 1 : 0) +
        species_habitat[, h] +
        (y == N_years ? species_habitat_change[, h] : rep_vector(0, N_species));

      vector[N_species] lambda_s = exp(log_expected);

      // --- Posterior predictive draws (for model checking and raw count summaries) ---
      for (s in 1:N_species)
        predicted_seed_counts[y, h, s] = neg_binomial_2_rng(lambda_s[s], phi[s]);

      // --- Expected species richness (analytical) ---
      // E[S] = sum_s P(Y_s > 0) = sum_s [1 - P(Y_s = 0)].
      // Derived from the NB pmf, so overdispersion (phi) is factored in:
      // a highly overdispersed species with a given lambda contributes less
      // to expected richness than a species with the same lambda but lower overdispersion.
      vector[N_species] p_present;
      for (s in 1:N_species)
        p_present[s] = 1.0 - exp(neg_binomial_2_lpmf(0 | lambda_s[s], phi[s]));
      richness[y, h] = sum(p_present);

      // --- Shannon entropy and Hill N1 (analytical) ---
      // H is computed from expected proportional abundances p_s = lambda_s / sum(lambda_s)
      // over all N_species. Using continuous lambda_s rather than integer predictive draws
      // avoids discretization artifacts (jagged posteriors) in H and Hill N1.
      // H is bounded above by log(N_species), so the evenness denominator uses
      // log(N_species) to guarantee J in [0, 1].
      real total_lambda = sum(lambda_s);
      real H = 0.0;
      for (s in 1:N_species) {
        real p = lambda_s[s] / total_lambda;
        if (p > 1e-15)
          H -= p * log(p);
      }

      // Hill N1: effective number of equally-abundant species on the Shannon scale.
      hill_N1[y, h] = exp(H);

      // Pielou's J: entropy relative to the maximum possible given the species pool.
      // Denominator is log(N_species) — the upper bound of H — ensuring J in [0, 1].
      // Interpretation: J = 1 when all species in the pool are equally abundant;
      // J = 0 when a single species dominates entirely.
      evenness[y, h] = (N_species > 1) ? H / log(N_species) : 0.0;
    }
  }
}
