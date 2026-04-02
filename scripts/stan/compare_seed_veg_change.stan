 // Hierarchical linked Gaussian models: vegetation and seed bank change
// Alec Chiono; alec.chiono@colorado.edu
// Version: b_veg[Kspecies, 2] (col1=mesic, col2=xeric) with correct habitat matching
//          + b_seed_mesic is now species-specific (vector[Kspecies])

data {
  // Veg data
  int<lower=1>               Kspecies;
  int<lower=1>               Nveg_plots;
  int<lower=1>               Nveg;
  array[Nveg] int<lower=1, upper=Kspecies> veg_species;
  array[Nveg] int<lower=1, upper=Nveg_plots> veg_plot;
  vector[Nveg]               veg_year;
  vector[Nveg]               veg_count;
  array[Nveg_plots] int<lower=0, upper=1> veg_habitat;

  // Seed data
  int<lower=1>               Nseed_plots;
  int<lower=1>               Nseed;
  vector[Nseed]              seed_change;
  array[Nseed] int<lower=1, upper=Kspecies> seed_species;
  array[Nseed] int<lower=1, upper=Nseed_plots> seed_plot;
  array[Nseed_plots] int<lower=0, upper=1> seed_habitat;
}

parameters {
  // Vegetation sub-model
  vector[Kspecies]             a_veg;
  matrix[Kspecies, 2]          b_veg;                  // col1=mesic, col2=xeric
  vector<lower=0>[Kspecies]    sigma_veg;
  matrix[Nveg_plots, Kspecies] z_veg_plot;
  vector<lower=0>[Kspecies]    sigma_veg_plot;

  // Seed sub-model
  vector[Kspecies]             a_seed;
  vector[Kspecies]             b_seed_mesic;            // ← now species-specific
  vector[Kspecies]             beta_seed_habitat_bveg;  // xeric shift on linkage slope
  vector[Kspecies]             beta_seed_habitat;       // xeric shift on intercept
  vector<lower=0>[Kspecies]    sigma_seed;
  vector[Nseed_plots]          z_seed_plot;
  real<lower=0>                sigma_seed_plot;
}

transformed parameters {
  // Veg plot random effects
  matrix[Nveg_plots, Kspecies] veg_plot_effect =
    z_veg_plot .* rep_matrix(sigma_veg_plot', Nveg_plots);

  vector[Nveg] veg_obs_plot_effect;
  for (i in 1:Nveg)
    veg_obs_plot_effect[i] = veg_plot_effect[veg_plot[i], veg_species[i]];

  // Habitat-specific veg slope per observation
  vector[Nveg] veg_slope;
  for (i in 1:Nveg)
    veg_slope[i] = b_veg[veg_species[i], veg_habitat[veg_plot[i]] + 1];

  // Seed plot random effects
  vector[Nseed_plots] seed_plot_effect = sigma_seed_plot * z_seed_plot;

  // Habitat-matched b_veg and habitat indicator per seed obs
  vector[Nseed] seed_b_veg;
  vector[Nseed] seed_hab_real;
  for (i in 1:Nseed) {
    int h = seed_habitat[seed_plot[i]];
    seed_b_veg[i]    = b_veg[seed_species[i], h + 1];
    seed_hab_real[i] = h * 1.0;
  }
}

model {
  // Priors — vegetation
  a_veg              ~ student_t(3, 0, 2.5);
  to_vector(b_veg)   ~ student_t(3, 0, 2.5);
  sigma_veg          ~ student_t(3, 0, 1);
  sigma_veg_plot     ~ student_t(3, 0, 1);
  to_vector(z_veg_plot) ~ std_normal();

  // Priors — seed
  a_seed                  ~ student_t(3, 0, 2.5);
  b_seed_mesic            ~ student_t(3, 0, 2.5);      // ← vector prior, same syntax
  beta_seed_habitat_bveg  ~ student_t(3, 0, 2.5);
  beta_seed_habitat       ~ student_t(3, 0, 2.5);
  sigma_seed              ~ student_t(3, 0, 1);
  sigma_seed_plot         ~ student_t(3, 0, 1);
  z_seed_plot             ~ std_normal();

  // Vegetation likelihood
  veg_count ~ normal(
    a_veg[veg_species] + veg_slope .* veg_year + veg_obs_plot_effect,
    sigma_veg[veg_species]
  );

  // Seed likelihood — b_seed_mesic now indexed by species
  seed_change ~ normal(
    a_seed[seed_species]
    + (b_seed_mesic[seed_species]                               // ← species-indexed
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
      + (b_seed_mesic[s]                                        // ← species-indexed
         + beta_seed_habitat_bveg[s] * h) * b_veg[s, h + 1]
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
