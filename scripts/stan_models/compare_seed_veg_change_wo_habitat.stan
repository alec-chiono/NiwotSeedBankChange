// Hierarchical linked Gaussian models: vegetation and seed bank change
// Alec Chiono; alec.chiono@colorado.edu
// Version: no habitat effects

data {
  // Veg data
  int<lower=1>               Kspecies;                     // number of species
  int<lower=1>               Nveg_plots;                   // number of veg plots
  int<lower=1>               Nveg;                         // total veg observations
  array[Nveg] int<lower=1, upper=Kspecies> veg_species;    // species ID per veg obs
  array[Nveg] int<lower=1, upper=Nveg_plots> veg_plot;     // plot ID per veg obs
  vector[Nveg]               veg_year;                     // year
  vector[Nveg]               veg_count;                    // veg response

  // Seed data
  int<lower=1>               Nseed_plots;                  // number of seed plots
  int<lower=1>               Nseed;                        // total seed observations
  vector[Nseed]              seed_change;                  // seed change response
  array[Nseed] int<lower=1, upper=Kspecies> seed_species;  // species ID per seed obs
  array[Nseed] int<lower=1, upper=Nseed_plots> seed_plot;  // plot ID per seed obs
}

parameters {
  // Vegetation sub-model
  vector[Kspecies]             a_veg;           // intercept per species
  vector[Kspecies]             b_veg;           // slope per species
  vector<lower=0>[Kspecies]    sigma_veg;       // species veg SD
  matrix[Nveg_plots, Kspecies] z_veg_plot;      // standardized veg plot RE per species
  vector<lower=0>[Kspecies]    sigma_veg_plot;  // SD of veg plot RE per species

  // Seed sub-model
  vector[Kspecies]           a_seed;            // intercept per species
  real                       b_seed;            // shared linkage slope
  vector<lower=0>[Kspecies]  sigma_seed;        // species seed SD
  vector[Nseed_plots]        z_seed_plot;       // standardized seed plot RE
  real<lower=0>              sigma_seed_plot;   // SD of seed plot RE
}

transformed parameters {
  // Veg plot random effects
  matrix[Nveg_plots, Kspecies] veg_plot_effect =
    z_veg_plot .* rep_matrix(sigma_veg_plot', Nveg_plots);

  vector[Nveg] veg_obs_plot_effect;
  for (i in 1:Nveg)
    veg_obs_plot_effect[i] = veg_plot_effect[veg_plot[i], veg_species[i]];

  // Seed plot random effects
  vector[Nseed_plots] seed_plot_effect = sigma_seed_plot * z_seed_plot;
}

model {
  // Priors
  a_veg ~ student_t(3, 0, 2.5);
  b_veg ~ student_t(3, 0, 2.5);
  sigma_veg ~ student_t(3, 0, 1);
  sigma_veg_plot ~ student_t(3, 0, 1);
  to_vector(z_veg_plot) ~ std_normal();

  sigma_seed_plot ~ student_t(3, 0, 1);
  z_seed_plot ~ std_normal();
  a_seed ~ student_t(3, 0, 2.5);
  b_seed ~ student_t(3, 0, 2.5);
  sigma_seed ~ student_t(3, 0, 1);

  // Vegetation likelihood
  veg_count ~ normal(
    a_veg[veg_species] + b_veg[veg_species] .* veg_year + veg_obs_plot_effect,
    sigma_veg[veg_species]
  );

  // Seed likelihood
  seed_change ~ normal(
    a_seed[seed_species] + b_seed * b_veg[seed_species] + seed_plot_effect[seed_plot],
    sigma_seed[seed_species]
  );
}

generated quantities {
  vector[Nseed] seed_change_rep;
  for (i in 1:Nseed) {
    int s = seed_species[i];
    int p = seed_plot[i];
    seed_change_rep[i] = normal_rng(
      a_seed[s] + b_seed * b_veg[s] + seed_plot_effect[p],
      sigma_seed[s]
    );
  }

  vector[Nveg] veg_count_rep;
  for (i in 1:Nveg) {
    veg_count_rep[i] = normal_rng(
      a_veg[veg_species[i]] + b_veg[veg_species[i]] * veg_year[i] + veg_obs_plot_effect[i],
      sigma_veg[veg_species[i]]
    );
  }
}
