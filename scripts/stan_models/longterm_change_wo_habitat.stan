// =============================================================================
// longterm_change_wo_habitat.stan
// =============================================================================
// Hierarchical linked Gaussian models: vegetation and seed bank change.
// Habitat-blind version: no habitat covariates in either submodel.
//
// PURPOSE
//   Jointly model (1) long-term vegetation cover change and (2) long-term seed
//   bank abundance change, linking the two via species-specific vegetation
//   temporal slopes (b_veg). This version pools estimation across habitats,
//   yielding a single community-level linkage coefficient per species rather
//   than habitat-stratified estimates. Intended for use as a null / reduced
//   model against longterm_change.stan.
//
// SUBMODEL 1 — VEGETATION
//   Response: veg_count (continuous; e.g. log-transformed cover)
//   Linear model per observation:
//     veg_count[i] ~ Normal(a_veg[s] + b_veg[s] * veg_year[i]
//                           + veg_plot_effect[p, s],
//                           sigma_veg[s])
//   where s = species, p = plot. Species-level intercepts and slopes are
//   drawn from hierarchical normal priors; plot effects are species-specific.
//
// SUBMODEL 2 — SEED BANK
//   Response: seed_change (continuous; change in log seed density, 1989→2023)
//   Linear model per observation:
//     seed_change[i] ~ Normal(a_seed[s] + b_seed[s] * b_veg[s]
//                             + seed_plot_effect[p],
//                             sigma_seed[s])
//   The key linkage term is b_veg[s]: the species-specific vegetation temporal
//   slope acts as a predictor of seed bank change. b_seed[s] is the
//   species-level coefficient scaling that linkage, drawn from a hierarchical
//   normal prior. Seed bank plot effects share a common SD across species.
//
// DIFFERENCES FROM longterm_change.stan
//   - No habitat arrays (veg_habitat, seed_habitat) in the data block.
//   - b_veg is a vector [Kspecies] rather than a matrix [Kspecies, 2];
//     a single temporal slope is estimated per species pooled across habitats.
//   - b_seed is species-specific (hierarchical) rather than a shared scalar;
//     the habitat-linkage interaction and direct habitat effect are absent.
//
// PARAMETERIZATION
//   All species-level random effects use non-centered (Matt trick)
//   parameterization. Residual and plot SDs use log-scale non-centering
//   (log(sigma) = log_mu + sigma_log * raw) to enforce positivity and improve
//   geometry when SDs might approach zero or vary substantially across species.
//
// AUTHOR
//   Alec Chiono; alec.chiono@colorado.edu
// =============================================================================


data {
  // ── Vegetation data ────────────────────────────────────────────────────────
  int<lower=1>                          Kspecies;   // number of species (shared across submodels)
  int<lower=1>                          Nveg_plots; // number of vegetation plots
  int<lower=1>                          Nveg;       // number of vegetation observations
  array[Nveg]      int<lower=1, upper=Kspecies>   veg_species; // species ID per veg observation
  array[Nveg]      int<lower=1, upper=Nveg_plots> veg_plot;    // plot ID per veg observation
  vector[Nveg]                          veg_year;   // year covariate (continuous or scaled)
  vector[Nveg]                          veg_count;  // vegetation response (continuous)

  // ── Seed bank data ─────────────────────────────────────────────────────────
  int<lower=1>                          Nseed_plots; // number of seed bank plots
  int<lower=1>                          Nseed;       // number of seed bank observations
  vector[Nseed]                         seed_change; // seed bank change response (continuous)
  array[Nseed]     int<lower=1, upper=Kspecies>    seed_species; // species ID per seed observation
  array[Nseed]     int<lower=1, upper=Nseed_plots> seed_plot;    // plot ID per seed observation
}


parameters {
  // ── Vegetation: species-level intercept ───────────────────────────────────
  // a_veg[s] ~ Normal(mu_a_veg, sigma_a_veg); non-centered below.
  real                  mu_a_veg;
  real<lower=0>         sigma_a_veg;
  vector[Kspecies]      a_veg_raw;       // a_veg = mu_a_veg + sigma_a_veg * a_veg_raw

  // ── Vegetation: species-level temporal slope ───────────────────────────────
  // A single slope per species, pooled across habitats.
  // b_veg[s] ~ Normal(mu_b_veg, sigma_b_veg); non-centered below.
  real                  mu_b_veg;
  real<lower=0>         sigma_b_veg;
  vector[Kspecies]      b_veg_raw;       // b_veg = mu_b_veg + sigma_b_veg * b_veg_raw

  // ── Vegetation: residual SD per species (log-scale non-centering) ─────────
  // log(sigma_veg[s]) = log_mu_sigma_veg + sigma_log_sigma_veg * sigma_veg_raw[s]
  // Log-scale non-centering enforces positivity and improves sampler geometry
  // when species-level SDs vary substantially.
  real                  log_mu_sigma_veg;
  real<lower=0>         sigma_log_sigma_veg;
  vector[Kspecies]      sigma_veg_raw;

  // ── Vegetation: plot SD per species (log-scale non-centering) ─────────────
  // Species-specific plot SDs allow some species to vary more across plots.
  real                  log_mu_sigma_veg_plot;
  real<lower=0>         sigma_log_sigma_veg_plot;
  vector[Kspecies]      sigma_veg_plot_raw;

  // ── Vegetation: plot random effects ───────────────────────────────────────
  // z_veg_plot[p, s] ~ std_normal(); scaled in transformed parameters by
  // sigma_veg_plot[s] to give the actual plot effect for species s at plot p.
  matrix[Nveg_plots, Kspecies] z_veg_plot;

  // ── Seed bank: species-level intercept ────────────────────────────────────
  real                  mu_a_seed;
  real<lower=0>         sigma_a_seed;
  vector[Kspecies]      a_seed_raw;

  // ── Seed bank: vegetation-linkage slope, species-specific ─────────────────
  // b_seed[s] scales how strongly species s's vegetation slope predicts its
  // seed bank change, pooled across habitats (no habitat interaction here).
  // b_seed[s] ~ Normal(mu_b_seed, sigma_b_seed); non-centered below.
  real                  mu_b_seed;
  real<lower=0>         sigma_b_seed;
  vector[Kspecies]      b_seed_raw;

  // ── Seed bank: residual SD per species (log-scale non-centering) ──────────
  real                  log_mu_sigma_seed;
  real<lower=0>         sigma_log_sigma_seed;
  vector[Kspecies]      sigma_seed_raw;

  // ── Seed bank: plot random effects ────────────────────────────────────────
  // Single shared sigma_seed_plot (not species-specific) because seed bank
  // plots may have fewer observations per species than vegetation plots,
  // making species-specific plot SDs poorly identified.
  vector[Nseed_plots]   z_seed_plot;
  real<lower=0>         sigma_seed_plot;
}


transformed parameters {
  // ── Recover centered vegetation species vectors ────────────────────────────
  vector[Kspecies] a_veg = mu_a_veg + sigma_a_veg * a_veg_raw;
  vector[Kspecies] b_veg = mu_b_veg + sigma_b_veg * b_veg_raw;

  // Residual and plot SDs on the natural scale via log-scale non-centering
  vector<lower=0>[Kspecies] sigma_veg =
    exp(log_mu_sigma_veg      + sigma_log_sigma_veg      * sigma_veg_raw);
  vector<lower=0>[Kspecies] sigma_veg_plot =
    exp(log_mu_sigma_veg_plot + sigma_log_sigma_veg_plot * sigma_veg_plot_raw);

  // ── Recover centered seed bank species vectors ─────────────────────────────
  vector[Kspecies] a_seed = mu_a_seed + sigma_a_seed * a_seed_raw;
  vector[Kspecies] b_seed = mu_b_seed + sigma_b_seed * b_seed_raw;

  vector<lower=0>[Kspecies] sigma_seed =
    exp(log_mu_sigma_seed + sigma_log_sigma_seed * sigma_seed_raw);

  // ── Plot random effects (vegetation) ──────────────────────────────────────
  // Element-wise scaling: each plot-species combination gets its own magnitude
  // governed by sigma_veg_plot[s], implemented as broadcasting via rep_matrix.
  matrix[Nveg_plots, Kspecies] veg_plot_effect =
    z_veg_plot .* rep_matrix(sigma_veg_plot', Nveg_plots);

  // Index plot effects to observations (avoids repeated indexing in model block)
  vector[Nveg] veg_obs_plot_effect;
  for (i in 1:Nveg)
    veg_obs_plot_effect[i] = veg_plot_effect[veg_plot[i], veg_species[i]];

  // Seed bank plot effects (scalar sigma; not species-specific)
  vector[Nseed_plots] seed_plot_effect = sigma_seed_plot * z_seed_plot;
}


model {
  // ── Hyperpriors ───────────────────────────────────────────────────────────
  // Weakly informative Student-t(3) for location hyperparameters; half-t for
  // scale hyperparameters. Log-scale priors for variance hyperparameters
  // correspond to a weakly informative prior on the natural-scale SDs,
  // with E[sigma] ≈ exp(-1) ≈ 0.37 under the Normal(-1, 1) prior.

  // Vegetation intercept hierarchy
  mu_a_veg                 ~ student_t(3, 0, 2.5);
  sigma_a_veg              ~ student_t(3, 0, 1);

  // Vegetation slope hierarchy
  mu_b_veg                 ~ student_t(3, 0, 2.5);
  sigma_b_veg              ~ student_t(3, 0, 1);

  // Vegetation residual SD hierarchy
  log_mu_sigma_veg         ~ normal(-1, 1);
  sigma_log_sigma_veg      ~ normal(0, 0.5);  // modest variation in residual SD across species

  // Vegetation plot SD hierarchy
  log_mu_sigma_veg_plot    ~ normal(-1, 1);
  sigma_log_sigma_veg_plot ~ normal(0, 0.5);

  // Seed bank intercept hierarchy
  mu_a_seed                ~ student_t(3, 0, 2.5);
  sigma_a_seed             ~ student_t(3, 0, 1);

  // Seed bank linkage slope hierarchy
  mu_b_seed                ~ student_t(3, 0, 2.5);
  sigma_b_seed             ~ student_t(3, 0, 1);

  // Seed bank residual SD hierarchy
  log_mu_sigma_seed        ~ normal(-1, 1);
  sigma_log_sigma_seed     ~ normal(0, 0.5);

  // Seed bank plot SD
  sigma_seed_plot          ~ student_t(3, 0, 1);

  // ── Non-centered raw parameters — all std_normal by construction ──────────
  a_veg_raw          ~ std_normal();
  b_veg_raw          ~ std_normal();
  sigma_veg_raw      ~ std_normal();
  sigma_veg_plot_raw ~ std_normal();
  a_seed_raw         ~ std_normal();
  b_seed_raw         ~ std_normal();
  sigma_seed_raw     ~ std_normal();

  // ── Plot random effect raw terms ──────────────────────────────────────────
  to_vector(z_veg_plot) ~ std_normal();
  z_seed_plot           ~ std_normal();

  // ── Likelihoods ───────────────────────────────────────────────────────────

  // Vegetation: species-specific residual SD, vectorized across observations
  veg_count ~ normal(
    a_veg[veg_species] + b_veg[veg_species] .* veg_year + veg_obs_plot_effect,
    sigma_veg[veg_species]
  );

  // Seed bank: species-level linkage, pooled across habitats
  //   mu[i] = a_seed[s] + b_seed[s] * b_veg[s] + seed_plot_effect[p]
  // b_seed[s] scales how strongly the vegetation temporal trend for species s
  // predicts its seed bank change; there is no habitat modulation of this
  // relationship in this version.
  seed_change ~ normal(
    a_seed[seed_species] + b_seed[seed_species] .* b_veg[seed_species]
    + seed_plot_effect[seed_plot],
    sigma_seed[seed_species]
  );
}


generated quantities {
  // ── Posterior predictive draws for retrodictive checks ────────────────────

  // Seed bank: replicated observations on the seed_change scale
  vector[Nseed] seed_change_rep;
  for (i in 1:Nseed) {
    int s = seed_species[i];
    int p = seed_plot[i];
    seed_change_rep[i] = normal_rng(
      a_seed[s] + b_seed[s] * b_veg[s] + seed_plot_effect[p],
      sigma_seed[s]
    );
  }

  // Vegetation: replicated observations on the veg_count scale
  vector[Nveg] veg_count_rep;
  for (i in 1:Nveg) {
    veg_count_rep[i] = normal_rng(
      a_veg[veg_species[i]] + b_veg[veg_species[i]] * veg_year[i] + veg_obs_plot_effect[i],
      sigma_veg[veg_species[i]]
    );
  }
}
