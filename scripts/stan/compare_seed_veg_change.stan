// Hierarchical linked Gaussian models: vegetation and seed bank change
// Alec Chiono; alec.chiono@colorado.edu
// Version: b_veg[Kspecies, 2] (col1=mesic, col2=xeric) with correct habitat matching

data {
  // Veg data
  int<lower=1>               Kspecies;                     // number of species
  int<lower=1>               Nveg_plots;                   // number of veg plots
  int<lower=1>               Nveg;                         // total veg observations
  array[Nveg] int<lower=1, upper=Kspecies> veg_species;    // species ID per veg obs
  array[Nveg] int<lower=1, upper=Nveg_plots> veg_plot;     // plot ID per veg obs
  vector[Nveg]               veg_year;                     // year (can be centered/scaled)
  vector[Nveg]               veg_count;                    // veg response
  array[Nveg_plots] int<lower=0, upper=1> veg_habitat;     // 0 = mesic, 1 = xeric (per veg plot)

  // Seed data
  int<lower=1>               Nseed_plots;                  // number of seed plots
  int<lower=1>               Nseed;                        // total seed observations
  vector[Nseed]              seed_change;                  // seed change response
  array[Nseed] int<lower=1, upper=Kspecies> seed_species;  // species ID per seed obs
  array[Nseed] int<lower=1, upper=Nseed_plots> seed_plot;  // plot ID per seed obs
  array[Nseed_plots] int<lower=0, upper=1> seed_habitat;   // 0 = mesic, 1 = xeric (per seed plot)
}

parameters {
  // Vegetation sub-model
  vector[Kspecies]  a_veg;                        // intercept per species
  matrix[Kspecies, 2]        b_veg;                        // slope per species & habitat (col1=mesic, col2=xeric)
  vector<lower=0>[Kspecies]  sigma_veg;                    // species veg SD
  matrix[Nveg_plots, Kspecies] z_veg_plot;                 // standardized veg plot RE per species
  vector<lower=0>[Kspecies]  sigma_veg_plot;               // SD of veg plot RE per species

  // Seed sub-model
  vector[Kspecies]           a_seed;                       // intercept per species
  real                       b_seed_mesic;                 // baseline linkage slope (mesic)
  vector[Kspecies]           beta_seed_habitat_bveg;       // habitat effect on linkage slope (xeric shift)
  vector[Kspecies]           beta_seed_habitat;            // habitat effect on seed intercept (xeric shift)
  vector<lower=0>[Kspecies]  sigma_seed;                   // species seed SD
  vector[Nseed_plots]        z_seed_plot;                  // standardized seed plot RE
  real<lower=0>              sigma_seed_plot;              // SD of seed plot RE
}

transformed parameters {
  // Veg plot random effects
  matrix[Nveg_plots, Kspecies] veg_plot_effect =
    z_veg_plot .* rep_matrix(sigma_veg_plot', Nveg_plots);

  vector[Nveg] veg_obs_plot_effect;
  for (i in 1:Nveg)
    veg_obs_plot_effect[i] = veg_plot_effect[veg_plot[i], veg_species[i]];

  // Veg slopes per observation: pick correct habitat-specific slope
  // habitat 0 -> column 1 (mesic), habitat 1 -> column 2 (xeric)
  vector[Nveg] veg_slope;
  for (i in 1:Nveg)
    veg_slope[i] = b_veg[veg_species[i], veg_habitat[veg_plot[i]] + 1];

  // Seed plot random effects
  vector[Nseed_plots] seed_plot_effect = sigma_seed_plot * z_seed_plot;

  // Precompute per-seed-obs quantities using habitat lookup
  vector[Nseed] seed_b_veg;         // habitat-matched b_veg for each seed obs
  vector[Nseed] seed_hab_real;       // habitat as real (0.0 or 1.0) for each seed obs
  for (i in 1:Nseed) {
    int h = seed_habitat[seed_plot[i]];
    seed_b_veg[i]    = b_veg[seed_species[i], h + 1];
    seed_hab_real[i] = h * 1.0;
  }
}

model {
  // Priors
  a_veg ~ student_t(3, 0, 2.5);
  to_vector(b_veg) ~ student_t(3, 0, 2.5);
  sigma_veg ~ student_t(3, 0, 1);
  sigma_veg_plot ~ student_t(3, 0, 1);
  to_vector(z_veg_plot) ~ std_normal();

  sigma_seed_plot ~ student_t(3, 0, 1);
  z_seed_plot ~ std_normal();
  a_seed ~ student_t(3, 0, 2.5);
  b_seed_mesic ~ student_t(3, 0, 2.5);
  beta_seed_habitat_bveg ~ student_t(3, 0, 2.5);
  beta_seed_habitat ~ student_t(3, 0, 2.5);
  sigma_seed ~ student_t(3, 0, 1);

  // Vegetation likelihood (habitat-specific slopes)
  veg_count ~ normal(
    a_veg[veg_species] + veg_slope .* veg_year + veg_obs_plot_effect,
    sigma_veg[veg_species]
  );

  // Seed likelihood: uses b_veg for the correct habitat
  // seed_habitat[seed_plot] is 0 or 1; +1 maps to b_veg col1 or col2
  seed_change ~ normal(
    a_seed[seed_species]
    + (b_seed_mesic
       + beta_seed_habitat_bveg[seed_species] .* seed_hab_real)
      .* seed_b_veg
    + seed_plot_effect[seed_plot]
    + beta_seed_habitat[seed_species] .* seed_hab_real,
    sigma_seed[seed_species]
  );
}

generated quantities {
  vector[Nseed] seed_change_rep;
  for (i in 1:Nseed) {
    int s = seed_species[i];
    int p = seed_plot[i];
    int h = seed_habitat[p];
    seed_change_rep[i] = normal_rng(
      a_seed[s]
      + (b_seed_mesic + beta_seed_habitat_bveg[s] * h) * b_veg[s, h + 1]
      + seed_plot_effect[p]
      + beta_seed_habitat[s] * h * 1.0,
      sigma_seed[s]
    );
  }

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

