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
  // Store log lambda and log likelihood for easier PPC later
  vector[N_obs] log_lambda_gq = log_lambda;
  vector[N_obs] log_lik;
  for (n in 1:N_obs) log_lik[n] = neg_binomial_2_log_lpmf(count[n] | log_lambda[n], phi);


  // Biodiversity indices by year and habitat
  array[N_years, N_habitats] int richness;     // count of predicted species above threshold
  array[N_years, N_habitats] real evenness;    // Pielou’s J' (H'/Hmax)
  array[N_years, N_habitats] real hill_N1;     // Hill number (exp(Shannon H'))
  // Predicted seed counts by year, habitat, and species (posterior draws)
  array[N_years, N_habitats, N_species] int predicted_seed_counts;

  for (y in 1:N_years) {
    for (h in 1:N_habitats) {
      // Expected log-abundance per species for this year–habitat combination
      vector[N_species] log_expected =
        mu_global +
        year_effect * (y == N_years ? 1 : 0) +
        habitat_effect * (h == 2 ? 1 : 0) +
        year_hab_int * (y == N_years ? 1 : 0) * (h == 2 ? 1 : 0) +
        species_habitat[, h] +
        (y == N_years ? species_habitat_change[, h] : rep_vector(0, N_species));

      vector[N_species] expected_abund = exp(log_expected);   // back-transform to abundance scale
      real total_abundance = sum(expected_abund);

      // Generate predicted seed counts from neg_binomial_2 (mean=exp(log_expected), phi)
      for (s in 1:N_species) {
        predicted_seed_counts[y, h, s] = neg_binomial_2_rng(exp(log_expected[s]), phi);
      }

      if (total_abundance > 0) {
        vector[N_species] p = expected_abund / total_abundance;  // relative abundance proportions

        // Species Richness
        richness[y, h] = 0;
        for (s in 1:N_species)
          if (expected_abund[s] >= .25)
            richness[y, h] += 1;

        // Compute Shannon’s H', which is used to calculate evenness and diversity
        vector[N_species] p_present = rep_vector(0.0, N_species);
        real total_present = 0.0;

        // Restrict to species judged “present” (> 0.5 abundance threshold)
        for (s in 1:N_species)
          if (expected_abund[s] > .25) {
            p_present[s] = p[s];
            total_present += p[s];
          }

        // Compute Shannon entropy with present species
        real H_present = 0.0;
        if (total_present > 0) {
          p_present = p_present / total_present;
          for (s in 1:N_species)
            if (p_present[s] > 1e-12)  // make sure don't take log of zero
              H_present += -p_present[s] * log(p_present[s]);
        }


        // Hill number N1 = exp(H’)
        hill_N1[y, h] = exp(H_present);

        // Evenness J’ = H’ / log(S)
        real S = richness[y, h];
        real H_max = (S > 1) ? log(S) : 0.0;
        evenness[y, h] = (H_max > 0) ? H_present / H_max : 0.0;

      } else { // If no abundance predicted, set biodiversity indices to zero
        richness[y, h] = 0;
        evenness[y, h] = 0.0;
        hill_N1[y, h] = 0.0;
      }
    }
  }
}
