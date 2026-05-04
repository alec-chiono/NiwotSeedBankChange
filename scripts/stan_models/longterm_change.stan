// =============================================================================
// longterm_change.stan
// =============================================================================
// Hierarchical linked Gaussian models: vegetation change and seed bank change.
//
// PURPOSE
//   Jointly model (1) long-term vegetation cover change and (2) long-term seed
//   bank abundance change, linking the two via species-specific vegetation
//   temporal slopes (b_veg). The seed bank submodel asks: does the direction
//   and magnitude of vegetation change over time predict seed bank change, and
//   does that linkage differ between xeric and mesic habitats?
//
// SUBMODEL 1 — VEGETATION
//   Response: veg_count (log-transformed cover or similar continuous measure)
//   Linear model per observation:
//     veg_count[i] ~ Normal(a_veg[s] + b_veg[s, h] * veg_year[i]
//                           + veg_plot_effect[p, s],
//                           sigma_veg[s])
//   where s = species, h = habitat (1 = xeric, 2 = mesic), p = plot.
//   Species-specific intercepts (a_veg) and habitat-specific slopes (b_veg)
//   both drawn from hierarchical normal priors. Plot random effects are
//   species-specific (matrix z_veg_plot) to allow plots to rank differently
//   across species.
//
// SUBMODEL 2 — SEED BANK
//   Response: seed_change (continuous change in log seed density, 1989→2023)
//   Linear model per observation:
//     seed_change[i] ~ Normal(
//       a_seed[s]
//       + (b_seed_mesic[s] + beta_seed_habitat_bveg[s] * habitat[p]) * b_veg[s, h]
//       + seed_plot_effect[p]
//       + beta_seed_habitat[s] * habitat[p],
//       sigma_seed[s])
//   The key linkage term is b_veg[s, h]: the species-specific vegetation
//   temporal slope acts as a predictor of seed bank change. The coefficients
//   b_seed_mesic and beta_seed_habitat_bveg allow this linkage to differ
//   between habitats (xeric reference, mesic additive shift).
//   beta_seed_habitat captures any habitat difference in seed change not
//   mediated through b_veg.
//
// HABITAT CODING
//   habitat is a binary integer: 0 = xeric (reference), 1 = mesic.
//   In the vegetation submodel, habitat + 1 indexes the b_veg column
//   (column 1 = xeric, column 2 = mesic).
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
  int<lower=1>                          Kspecies;    // number of species (shared across submodels)
  int<lower=1>                          Nveg_plots;  // number of vegetation plots
  int<lower=1>                          Nveg;        // number of vegetation observations
  array[Nveg]      int<lower=1, upper=Kspecies>   veg_species;  // species ID per veg observation
  array[Nveg]      int<lower=1, upper=Nveg_plots> veg_plot;     // plot ID per veg observation
  vector[Nveg]                          veg_year;    // year covariate (continuous or scaled)
  vector[Nveg]                          veg_count;   // vegetation response (continuous)
  array[Nveg_plots] int<lower=0, upper=1>          veg_habitat; // habitat per veg plot (0=xeric, 1=mesic)

  // ── Seed bank data ─────────────────────────────────────────────────────────
  int<lower=1>                          Nseed_plots; // number of seed bank plots
  int<lower=1>                          Nseed;       // number of seed bank observations
  vector[Nseed]                         seed_change; // seed bank change response (continuous)
  array[Nseed]     int<lower=1, upper=Kspecies>    seed_species; // species ID per seed observation
  array[Nseed]     int<lower=1, upper=Nseed_plots> seed_plot;    // plot ID per seed observation
  array[Nseed_plots] int<lower=0, upper=1>         seed_habitat; // habitat per seed plot (0=xeric, 1=mesic)
}


parameters {
  // ── Vegetation: species-level intercept ───────────────────────────────────
  // a_veg[s] ~ Normal(mu_a_veg, sigma_a_veg); non-centered below.
  real                  mu_a_veg;
  real<lower=0>         sigma_a_veg;
  vector[Kspecies]      a_veg_raw;        // a_veg = mu_a_veg + sigma_a_veg * a_veg_raw

  // ── Vegetation: species × habitat temporal slopes ─────────────────────────
  // b_veg[s, h] gives the species-specific vegetation trend for each habitat.
  // Two habitat columns allow xeric and mesic trajectories to be estimated
  // independently for each species rather than fitting a single shared slope
  // modified by an additive interaction.
  vector[2]             mu_b_veg;         // mean slope per habitat
  vector<lower=0>[2]    sigma_b_veg;      // SD of slopes per habitat
  matrix[Kspecies, 2]   b_veg_raw;        // b_veg[,h] = mu_b_veg[h] + sigma_b_veg[h] * b_veg_raw[,h]

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

  // ── Seed bank: vegetation-linkage slope (mesic reference) ─────────────────
  // b_seed_mesic[s] scales how strongly species s's vegetation slope predicts
  // its seed bank change in mesic habitat (the reference level for the linkage).
  real                  mu_b_seed_mesic;
  real<lower=0>         sigma_b_seed_mesic;
  vector[Kspecies]      b_seed_mesic_raw;

  // ── Seed bank: habitat shift on vegetation-linkage slope ──────────────────
  // beta_seed_habitat_bveg[s] is the additive xeric–mesic difference in how
  // strongly b_veg predicts seed bank change for species s.
  real                  mu_beta_seed_habitat_bveg;
  real<lower=0>         sigma_beta_seed_habitat_bveg;
  vector[Kspecies]      beta_seed_habitat_bveg_raw;

  // ── Seed bank: direct habitat effect on seed change ───────────────────────
  // beta_seed_habitat[s] captures any habitat difference in seed bank change
  // that is not mediated through b_veg (i.e., not explained by vegetation trend).
  real                  mu_beta_seed_habitat;
  real<lower=0>         sigma_beta_seed_habitat;
  vector[Kspecies]      beta_seed_habitat_raw;

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

  // Habitat-specific slopes recovered column-wise
  matrix[Kspecies, 2] b_veg;
  for (h in 1:2)
    b_veg[, h] = mu_b_veg[h] + sigma_b_veg[h] * b_veg_raw[, h];

  // Residual and plot SDs on natural scale via log-scale non-centering
  vector<lower=0>[Kspecies] sigma_veg =
    exp(log_mu_sigma_veg      + sigma_log_sigma_veg      * sigma_veg_raw);
  vector<lower=0>[Kspecies] sigma_veg_plot =
    exp(log_mu_sigma_veg_plot + sigma_log_sigma_veg_plot * sigma_veg_plot_raw);

  // ── Recover centered seed bank species vectors ────────────────────────────
  vector[Kspecies] a_seed =
    mu_a_seed + sigma_a_seed * a_seed_raw;
  vector[Kspecies] b_seed_mesic =
    mu_b_seed_mesic + sigma_b_seed_mesic * b_seed_mesic_raw;
  vector[Kspecies] beta_seed_habitat_bveg =
    mu_beta_seed_habitat_bveg + sigma_beta_seed_habitat_bveg * beta_seed_habitat_bveg_raw;
  vector[Kspecies] beta_seed_habitat =
    mu_beta_seed_habitat + sigma_beta_seed_habitat * beta_seed_habitat_raw;

  vector<lower=0>[Kspecies] sigma_seed =
    exp(log_mu_sigma_seed + sigma_log_sigma_seed * sigma_seed_raw);

  // ── Plot random effects (vegetation) ─────────────────────────────────────
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

  // ── Habitat-matched vegetation slope per observation ─────────────────────
  // Select the appropriate habitat column (1=xeric, 2=mesic) of b_veg for
  // each vegetation observation based on that observation's plot habitat.
  vector[Nveg] veg_slope;
  for (i in 1:Nveg)
    veg_slope[i] = b_veg[veg_species[i], veg_habitat[veg_plot[i]] + 1];

  // ── Habitat-matched b_veg and habitat indicator per seed observation ──────
  // Pre-computed here to keep the model block likelihood statement vectorized.
  // seed_b_veg[i]: vegetation slope for species s in habitat h (the linkage predictor)
  // seed_hab_real[i]: habitat as real (0.0 or 1.0) for arithmetic use in model block
  vector[Nseed] seed_b_veg;
  vector[Nseed] seed_hab_real;
  for (i in 1:Nseed) {
    int h = seed_habitat[seed_plot[i]];
    seed_b_veg[i]    = b_veg[seed_species[i], h + 1];
    seed_hab_real[i] = h * 1.0;
  }
}


model {
  // ── Hyperpriors ───────────────────────────────────────────────────────────
  // Weakly informative Student-t(3) for location hyperparameters; half-t for
  // scale hyperparameters. Log-scale priors for variance hyperparameters
  // correspond to a weakly informative prior on the natural scale SDs.

  // Vegetation intercept hierarchy
  mu_a_veg                 ~ student_t(3, 0, 2.5);
  sigma_a_veg              ~ student_t(3, 0, 1);

  // Vegetation slope hierarchy
  mu_b_veg                 ~ student_t(3, 0, 2.5);
  sigma_b_veg              ~ student_t(3, 0, 1);

  // Vegetation residual SD hierarchy: prior on log scale, E[sigma] ≈ exp(-1) ≈ 0.37
  log_mu_sigma_veg         ~ normal(-1, 1);
  sigma_log_sigma_veg      ~ normal(0, 0.5);  // modest variation in residual SD across species

  // Vegetation plot SD hierarchy
  log_mu_sigma_veg_plot    ~ normal(-1, 1);
  sigma_log_sigma_veg_plot ~ normal(0, 0.5);

  // Seed bank intercept hierarchy
  mu_a_seed                ~ student_t(3, 0, 2.5);
  sigma_a_seed             ~ student_t(3, 0, 1);

  // Seed bank linkage slope hierarchy (mesic reference)
  mu_b_seed_mesic          ~ student_t(3, 0, 2.5);
  sigma_b_seed_mesic       ~ student_t(3, 0, 1);

  // Seed bank habitat shift on linkage slope
  mu_beta_seed_habitat_bveg   ~ student_t(3, 0, 2.5);
  sigma_beta_seed_habitat_bveg ~ student_t(3, 0, 1);

  // Seed bank direct habitat effect
  mu_beta_seed_habitat     ~ student_t(3, 0, 2.5);
  sigma_beta_seed_habitat  ~ student_t(3, 0, 1);

  // Seed bank residual SD hierarchy
  log_mu_sigma_seed        ~ normal(-1, 1);
  sigma_log_sigma_seed     ~ normal(0, 0.5);

  // Seed bank plot SD
  sigma_seed_plot          ~ student_t(3, 0, 1);

  // ── Non-centered raw parameters — all std_normal by construction ──────────
  a_veg_raw                    ~ std_normal();
  to_vector(b_veg_raw)         ~ std_normal();
  sigma_veg_raw                ~ std_normal();
  sigma_veg_plot_raw           ~ std_normal();
  a_seed_raw                   ~ std_normal();
  b_seed_mesic_raw             ~ std_normal();
  beta_seed_habitat_bveg_raw   ~ std_normal();
  beta_seed_habitat_raw        ~ std_normal();
  sigma_seed_raw               ~ std_normal();

  // ── Plot random effect raw terms ──────────────────────────────────────────
  to_vector(z_veg_plot) ~ std_normal();
  z_seed_plot           ~ std_normal();

  // ── Likelihoods ───────────────────────────────────────────────────────────

  // Vegetation: species-specific residual SD, vectorized across observations
  veg_count ~ normal(
    a_veg[veg_species] + veg_slope .* veg_year + veg_obs_plot_effect,
    sigma_veg[veg_species]
  );

  // Seed bank: full linkage structure
  //   mu[i] = a_seed[s]
  //           + (b_seed_mesic[s] + beta_seed_habitat_bveg[s] * habitat[p]) * b_veg[s, h]
  //           + seed_plot_effect[p]
  //           + beta_seed_habitat[s] * habitat[p]
  // The first three terms model habitat-modulated vegetation-driven change;
  // the last term captures residual habitat differences in seed bank change.
  seed_change ~ normal(
    a_seed[seed_species]
    + (b_seed_mesic[seed_species]
       + beta_seed_habitat_bveg[seed_species] .* seed_hab_real)
      .* seed_b_veg
    + seed_plot_effect[seed_plot]
    + beta_seed_habitat[seed_species] .* seed_hab_real,
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
    int h = seed_habitat[p];
    seed_change_rep[i] = normal_rng(
      a_seed[s]
      + (b_seed_mesic[s]
         + beta_seed_habitat_bveg[s] * h) * b_veg[s, h + 1]
      + seed_plot_effect[p]
      + beta_seed_habitat[s] * h * 1.0,
      sigma_seed[s]
    );
  }

  // Vegetation: replicated observations on the veg_count scale
  vector[Nveg] veg_count_rep;
  for (i in 1:Nveg) {
    veg_count_rep[i] = normal_rng(
      a_veg[veg_species[i]]
      + veg_slope[i] * veg_year[i]
      + veg_obs_plot_effect[i],
      sigma_veg[veg_species[i]]
    );
  }
}
