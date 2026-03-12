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
  vector[N_years] year_effect_raw;             // standardized year effects
  vector[N_habitats] habitat_effect_raw;       // standardized habitat effects
  matrix[N_years, N_habitats] year_hab_int_raw;// standardized year-by-habitat interactions
  vector[N_plots] plot_effect_raw;             // standardized plot effects
  real<lower=0> sigma_year;                    // SD of year effects
  real<lower=0> sigma_habitat;                 // SD of habitat effects
  real<lower=0> sigma_interaction;             // SD of year:habitat interaction effects
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
  vector[N_years] year_effect = sigma_year * year_effect_raw;
  vector[N_habitats] habitat_effect = sigma_habitat * habitat_effect_raw;
  matrix[N_years, N_habitats] year_hab_int = sigma_interaction * year_hab_int_raw;
  vector[N_plots] plot_effect = sigma_plot * plot_effect_raw;
  matrix[N_species, N_habitats] species_habitat = sigma_species_habitat * species_habitat_raw;
  matrix[N_species, N_habitats] species_habitat_change = sigma_species_habitat_change * species_habitat_change_raw;


  // Construct the expected log-abundance for each observation
  for (i in 1:N_obs) {
    // Temporal change: apply only to the most recent year
    real temporal_change = (year_idx[i] == N_years) ? species_habitat_change[species[i], habitat_idx[i]] : 0;

    log_lambda[i] =
      mu_global +
      year_effect[year_idx[i]] +
      habitat_effect[habitat_idx[i]] +
      year_hab_int[year_idx[i], habitat_idx[i]] +
      plot_effect[plot[i]] +
      species_habitat[species[i], habitat_idx[i]] +
      temporal_change;
  }
}

model {
  // Priors (weakly informative)
  // Global intercept
  mu_global ~ student_t(3, 0, 2.5);
  // SDs slightly tighter
  sigma_year ~ normal(0, 0.5);
  sigma_habitat ~ normal(0, 2);
  sigma_interaction ~ normal(0, 0.5); // interactions even tighter
  sigma_plot ~ normal(0, 0.5);
  sigma_species_habitat ~ normal(0, 0.5);
  sigma_species_habitat_change ~ normal(0, 0.25); // interactions even tighter
  // All "raw" terms standard normal
  year_effect_raw ~ std_normal();
  habitat_effect_raw ~ std_normal();
  to_vector(year_hab_int_raw) ~ std_normal();
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
        year_effect[y] +
        habitat_effect[h] +
        year_hab_int[y, h] +
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
          if (expected_abund[s] >= 0.5)
            richness[y, h] += 1;

        // Compute Shannon’s H', which is used to calculate evenness and diversity
        vector[N_species] p_present = rep_vector(0.0, N_species);
        real total_present = 0.0;

        // Restrict to species judged “present” (> 0.5 abundance threshold)
        for (s in 1:N_species)
          if (expected_abund[s] > 0.5) {
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
