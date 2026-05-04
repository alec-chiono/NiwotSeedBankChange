// =============================================================================
// community_indices.stan
// =============================================================================
// Hierarchical negative binomial model for seed bank community composition.
//
// PURPOSE
//   Estimate species-level seed counts across a two-year (1989, 2023),
//   two-habitat (xeric, mesic) design, then derive posterior distributions
//   of three biodiversity indices per year × habitat combination:
//     - Expected species richness (E[S])
//     - Hill number N1 (Shannon effective species, exp(H))
//     - Pielou's evenness (J = H / log(S_pool))
//
// RESPONSE VARIABLE
//   Seed counts per species per plot (integer, overdispersed). Modeled with
//   a negative binomial likelihood (NB2 parameterization) to accommodate
//   zero inflation and variance > mean typical of seed bank data.
//
// FIXED EFFECTS (log scale)
//   - Global intercept
//   - Year effect (binary: 1989 = reference, 2023 = treatment)
//   - Habitat effect (binary: xeric = reference, mesic = treatment)
//   - Year × habitat interaction
//
// RANDOM EFFECTS
//   - Plot-level intercept (shared across species within a plot)
//   - Species × habitat baseline interaction (allows each species to differ
//     in abundance across habitats)
//   - Species × habitat temporal change (allows each species to respond
//     differently to the 1989→2023 transition within each habitat)
//
// OVERDISPERSION
//   Species-specific phi drawn from a log-normal hyperprior, allowing partial
//   pooling of dispersion estimates across species.
//
// PARAMETERIZATION
//   All random effects use non-centered (Matt trick) parameterization to
//   improve geometry and sampler efficiency in hierarchical models.
//
// GENERATED QUANTITIES
//   - Pointwise posterior predictive draws (obs-level, for LOO/retrodiction)
//   - Predicted seed counts per year × habitat × species (plot-marginalized)
//   - Biodiversity indices computed analytically from continuous lambda_s
//     rather than from integer predictive draws to avoid discretization
//     artifacts in H and Hill N1 (see GQ block comments).
//
// AUTHOR
//   Alec Chiono; alec.chiono@colorado.edu
// =============================================================================


data {
  int<lower=1> N_obs;                    // total number of observations
  int<lower=1> N_species;                // number of species
  int<lower=1> N_plots;                  // number of sampling plots
  int<lower=1> N_years;                  // number of years in the study
  int<lower=1> N_habitats;               // number of habitat types
  array[N_obs] int<lower=0> count;       // observed seed count (response variable)
  array[N_obs] int<lower=1> species;     // species ID for each observation
  array[N_obs] int<lower=1> plot;        // plot ID for each observation
  array[N_obs] int<lower=1> year_idx;    // year index (1 = 1989, N_years = 2023)
  array[N_obs] int<lower=1> habitat_idx; // habitat index (1 = xeric, 2 = mesic)
}


parameters {
  // ── Fixed effects ──────────────────────────────────────────────────────────
  real mu_global;      // global intercept (log scale)
  real year_effect;    // year contrast: 2023 vs. 1989
  real habitat_effect; // habitat contrast: mesic vs. xeric
  real year_hab_int;   // year × habitat interaction

  // ── Plot random effects (non-centered) ────────────────────────────────────
  // Actual plot effects recovered in transformed parameters as
  // plot_effect = sigma_plot * plot_effect_raw.
  vector[N_plots] plot_effect_raw;
  real<lower=0> sigma_plot;            // SD of plot effects

  // ── Species × habitat random effects (non-centered) ───────────────────────
  // Baseline species–habitat affinity and temporal change in that affinity.
  // Separate sigma hyperparameters allow the model to learn independently how
  // variable species are in their baseline habitat preferences vs. how variable
  // they are in their response to the 1989→2023 transition.
  matrix[N_species, N_habitats] species_habitat_raw;
  matrix[N_species, N_habitats] species_habitat_change_raw;
  real<lower=0> sigma_species_habitat;        // SD of baseline species–habitat effects
  real<lower=0> sigma_species_habitat_change; // SD of temporal change in species–habitat effects

  // ── Species-specific overdispersion (non-centered, log scale) ─────────────
  // phi_s ~ LogNormal(mu_log_phi, sigma_log_phi): partial pooling across species.
  // NB2 variance = lambda + lambda^2 / phi, so larger phi → less overdispersion.
  vector[N_species] phi_raw;
  real mu_log_phi;                // mean log-phi across species
  real<lower=0> sigma_log_phi;   // SD of log-phi across species
}


transformed parameters {
  // ── Recover centered random effects ───────────────────────────────────────
  vector[N_plots] plot_effect = sigma_plot * plot_effect_raw;
  matrix[N_species, N_habitats] species_habitat        = sigma_species_habitat        * species_habitat_raw;
  matrix[N_species, N_habitats] species_habitat_change = sigma_species_habitat_change * species_habitat_change_raw;

  // ── Species-specific overdispersion on the natural scale ──────────────────
  vector<lower=0>[N_species] phi = exp(mu_log_phi + sigma_log_phi * phi_raw);

  // ── Log-scale linear predictor ────────────────────────────────────────────
  // Temporal change applied only in the most recent year (year_idx == N_years).
  // Year and habitat contrasts use indicator arithmetic rather than separate
  // index arrays to keep the predictor construction compact and readable.
  vector[N_obs] log_lambda;
  for (i in 1:N_obs) {
    real temporal_change = (year_idx[i] == N_years)
                           ? species_habitat_change[species[i], habitat_idx[i]]
                           : 0;

    log_lambda[i] =
      mu_global +
      year_effect    * (year_idx[i] == N_years ? 1 : 0) +
      habitat_effect * (habitat_idx[i] == 2    ? 1 : 0) +
      year_hab_int   * (year_idx[i] == N_years ? 1 : 0) * (habitat_idx[i] == 2 ? 1 : 0) +
      plot_effect[plot[i]] +
      species_habitat[species[i], habitat_idx[i]] +
      temporal_change;
  }
}


model {
  // ── Priors ─────────────────────────────────────────────────────────────────
  // Weakly informative Student-t(3) priors on fixed effects. Heavier tails
  // than a normal allow occasional larger effects without strong regularization.

  // Fixed effects — moderate scale reflecting log-count data
  mu_global      ~ student_t(3, 0, 2.5);
  year_effect    ~ student_t(3, 0, 2.5);
  habitat_effect ~ student_t(3, 0, 2.5);
  year_hab_int   ~ student_t(3, 0, 1);   // tighter: interactions less likely to be large

  // SD hyperpriors — half-normal; tighter than fixed effects to encourage pooling
  sigma_plot                  ~ normal(0, 0.5);
  sigma_species_habitat       ~ normal(0, 0.5);
  sigma_species_habitat_change ~ normal(0, 0.25); // extra regularization on change terms

  // Non-centered raw parameters — standard normal by construction
  plot_effect_raw                    ~ std_normal();
  to_vector(species_habitat_raw)     ~ std_normal();
  to_vector(species_habitat_change_raw) ~ std_normal();

  // Overdispersion hyperpriors
  phi_raw       ~ std_normal();
  mu_log_phi    ~ normal(0, 1);    // weakly informative on log scale
  sigma_log_phi ~ normal(0, 0.5); // modest species-to-species variation in dispersion

  // ── Likelihood ────────────────────────────────────────────────────────────
  // NB2 parameterized via log_lambda for numerical stability;
  // phi is indexed by species rather than observation to enforce partial pooling.
  for (i in 1:N_obs)
    count[i] ~ neg_binomial_2_log(log_lambda[i], phi[species[i]]);
}


generated quantities {
  // ── Pointwise posterior predictive draws (obs-matched) ────────────────────
  // Used for LOO-CV and retrodictive checks; each draw matched 1:1 to an
  // observed count so residual patterns can be examined on the original scale.
  array[N_obs] int pred_count;
  for (i in 1:N_obs)
    pred_count[i] = neg_binomial_2_log_rng(log_lambda[i], phi[species[i]]);

  // ── Community-level predicted counts (year × habitat × species) ───────────
  // Plot effects are marginalized out: because plot effects are drawn from a
  // mean-zero normal prior, their expectation is zero and omitting them gives
  // the marginal (population-average) community prediction. These draws
  // represent the expected community composition at a "typical" plot.
  array[N_years, N_habitats, N_species] int predicted_seed_counts;

  // ── Biodiversity indices (year × habitat) ─────────────────────────────────
  array[N_years, N_habitats] real richness;  // E[S]: expected species richness
  array[N_years, N_habitats] real hill_N1;   // Hill N1: exp(H), Shannon effective species
  array[N_years, N_habitats] real evenness;  // Pielou's J = H / log(S_pool)

  for (y in 1:N_years) {
    for (h in 1:N_habitats) {

      // Expected log-abundance per species (plot effects marginalized out)
      vector[N_species] log_expected =
        mu_global +
        year_effect    * (y == N_years ? 1 : 0) +
        habitat_effect * (h == 2       ? 1 : 0) +
        year_hab_int   * (y == N_years ? 1 : 0) * (h == 2 ? 1 : 0) +
        species_habitat[, h] +
        (y == N_years ? species_habitat_change[, h] : rep_vector(0, N_species));

      vector[N_species] lambda_s = exp(log_expected);

      // Posterior predictive draws for raw count summaries and model checking
      for (s in 1:N_species)
        predicted_seed_counts[y, h, s] = neg_binomial_2_rng(lambda_s[s], phi[s]);

      // ── Expected species richness: E[S] = Σ_s P(Y_s > 0) ─────────────────
      // = Σ_s [1 − P(Y_s = 0 | lambda_s, phi_s)].
      // Derived analytically from the NB pmf so that overdispersion (phi) is
      // properly accounted for: a highly overdispersed species contributes less
      // to E[S] than an equally abundant species with lower overdispersion.
      vector[N_species] p_present;
      for (s in 1:N_species)
        p_present[s] = 1.0 - exp(neg_binomial_2_lpmf(0 | lambda_s[s], phi[s]));
      richness[y, h] = sum(p_present);

      // ── Shannon entropy and Hill N1 ───────────────────────────────────────
      // H is computed from the continuous expected proportional abundances
      // p_s = lambda_s / Σ lambda_s rather than from integer predictive draws.
      // Using continuous lambda avoids discretization artifacts (jagged
      // posterior distributions) in H and Hill N1 that arise when many species
      // have low predicted counts and their presence/absence fluctuates
      // stochastically across draws.
      //
      // The guard (p > 1e-15) prevents log(0) for species with negligible
      // expected abundance; in practice this matters only for species whose
      // lambda_s rounds to machine zero.
      real total_lambda = sum(lambda_s);
      real H = 0.0;
      for (s in 1:N_species) {
        real p = lambda_s[s] / total_lambda;
        if (p > 1e-15)
          H -= p * log(p);
      }

      // Hill N1: effective number of equally-abundant species on the Shannon scale
      hill_N1[y, h] = exp(H);

      // Pielou's J: H relative to its upper bound log(S_pool).
      // Using log(N_species) as denominator — the maximum H achievable given
      // the observed species pool — guarantees J ∈ [0, 1].
      // J = 1: all species equally abundant; J → 0: single species dominates.
      // Returns 0 for a single-species pool to avoid division by zero.
      evenness[y, h] = (N_species > 1) ? H / log(N_species) : 0.0;
    }
  }
}
