// Heirarchical linked gaussian models to relationship between vegetation and seed bank change
// Alec Chiono; alec.chiono@colorado.edu

data {
  // Veg data
  int<lower=1>               Kspecies;                    // number of species
  int<lower=1>               Nveg_plots;                  // number of plots
  int<lower=1>               Nveg;                        // total number of veg observations
  array[Nveg] int<lower=1, upper=Kspecies> veg_species;   // species ID for each veg observation
  array[Nveg] int<lower=1, upper=Nveg_plots> veg_plot;    // plot ID for each veg observation
  vector[Nveg]               veg_year;                    // year value for each veg observation
  vector[Nveg]               veg_count;                   // vegetation count for each veg observation
  // Seed Data
  int<lower=1>               Nseed_plots;                 // number of seed plots
  int<lower=1>               Nseed;                       // total number of seed observations
  vector[Nseed]              seed_change;                 // observed change in seed
  array[Nseed] int<lower=1, upper=Kspecies> seed_species; // species ID for seed observation
  array[Nseed] int<lower=1, upper=Nseed_plots> seed_plot; // plot ID for seed observation
}

parameters {
  // Vegetation sub-model parameters
  vector<lower=0>[Kspecies]  a_veg;     // intercept per species
  vector[Kspecies]           b_veg;     // slope per species
  vector<lower=0>[Kspecies]  sigma_veg; // SD per species
  // Veg plot random intercepts
  matrix[Nveg_plots, Kspecies]   z_veg_plot;      // standardized veg plot effects per species
  vector<lower=0>[Kspecies]      sigma_veg_plot;  // SD of veg plot effects per species

  // Seed sub-model parameters
  vector[Kspecies]           a_seed;        // intercept per species
  real                       b_seed;        // one shared slope
  vector<lower=0>[Kspecies]  sigma_seed;    // SD per species
  // Seed plot random intercepts
  vector[Nseed_plots]            z_seed_plot;     // standardized seed plot effects (shared)
  real<lower=0>                  sigma_seed_plot; // SD of seed plot effects (shared)
}

transformed parameters {
  // Actual plot effects for veg data
  matrix[Nveg_plots, Kspecies] veg_plot_effect = z_veg_plot .* rep_matrix(sigma_veg_plot', Nveg_plots);
  // Pre-compute veg plot effects for each observation
  vector[Nveg] veg_obs_plot_effect;
  for (i in 1:Nveg) veg_obs_plot_effect[i] = veg_plot_effect[veg_plot[i], veg_species[i]];

  // Actual plot effects for seed data
  vector[Nseed_plots] seed_plot_effect = sigma_seed_plot * z_seed_plot;
}

model {
  // Priors (weakly informative)
  a_veg ~ student_t(3, 0, 2.5);
  b_veg ~ student_t(3, 0, 2.5);
  sigma_veg ~ student_t(3, 0, 2.5);
  sigma_veg_plot ~ student_t(3, 0, 1);
  to_vector(z_veg_plot) ~ std_normal();
  sigma_seed_plot ~ student_t(3, 0, 1);
  z_seed_plot ~ std_normal();
  a_seed ~ student_t(3, 0, 2.5);
  b_seed ~ student_t(3, 0, 2.5);
  sigma_seed ~ student_t(3, 0, 2.5);

  // Likelihoods
  veg_count ~ normal(a_veg[veg_species] + b_veg[veg_species] .* veg_year + veg_obs_plot_effect, sigma_veg[veg_species]);
  seed_change ~ normal(a_seed[seed_species] + b_seed * b_veg[seed_species] + seed_plot_effect[seed_plot], sigma_seed[seed_species]);
}
