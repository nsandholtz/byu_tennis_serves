data {
  int<lower=0> N_obs;               // total number of uncensored observations
  int<lower=0> N_cens;              // total number of censored observations
  int<lower=0> N_p;                 // number of players
  int<lower=1> N_i;                 // number of serve periods
  int<lower=1> N_j;                 // Number of regions
  
  array[N_obs] row_vector[2] z;             // Observed (x, y) serve locations
  array[N_p, N_j] row_vector[2] m;          // Player/region-specific (x,y) target locations
  array[N_p, N_i, N_j] real a;              // player-serve period-region lower x-truncation max value

  array[N_obs] int<lower=1, upper=N_p> z_ind_p;  // Player ID for each obs serve
  array[N_obs] int<lower=1, upper=N_i> z_ind_i;  // Serve period for each obs
  array[N_obs] int<lower=1, upper=N_j> z_ind_j;  // Region ID for each obs
 
  array[N_cens] int<lower=1, upper=N_p> c_ind_p; // Player ID for censored
  array[N_cens] int<lower=1, upper=N_i> c_ind_i; // Serve period for censored
  array[N_cens] int<lower=1, upper=N_j> c_ind_j; // Region ID for censored
}

parameters {
  array[N_p, N_i, N_j] real rho_TRANS;                // transformed player-period-region correlations
  array[N_j] real rho_g_TRANS;                        // transformed global region correlations
  array[N_p, N_i, N_j] vector[2] tau_TRANS;           // transformed player-period-region scales
  array[N_p] vector<lower=0>[2] tau_g;                // player level scales

  array[N_p, N_i, N_j] vector[2] mu;                   // Mean vector per player-period-region
  array[N_p, N_i, N_j] real<lower=0, upper=a> alpha;   // Latent truncation threshold per player-period-region
  array[N_i] real<lower=0> alpha_g;                    // Latent truncation per serve period

  real<lower=0> sig_mu;                                // SD Hyperprior on mu
  real<lower=0> sig_alpha;                             // SD Hyperprior on alpha
  real<lower=0> sig_rho_TRANS;                         // SD Hyperprior on rho
  real<lower=0> sig_tau_TRANS;                         // SD Hyperprior on tau
}

transformed parameters {
  array[N_p, N_i, N_j] vector<lower=0>[2] tau;  
  array[N_p, N_i, N_j] real rho;
  array[N_j] real rho_g;
  array[N_p, N_i, N_j] corr_matrix[2] Omega;
  
  // recover rho_g
  for (j in 1:N_j) rho_g[j] = tanh(rho_g_TRANS[j]);
  
  // create correlation matrix
  for (p in 1:N_p) {
    for (i in 1:N_i) {
      for (j in 1:N_j) {
        rho[p, i, j] = tanh(rho_g_TRANS[j] + sig_rho_TRANS * rho_TRANS[p, i, j]);
        Omega[p, i, j][1, 1] = 1;
        Omega[p, i, j][2, 2] = 1;
        Omega[p, i, j][1, 2] = rho[p, i, j];
        Omega[p, i, j][2, 1] = rho[p, i, j];
      }
    }
  }
  
  // transform log-scale non-centered tau
  for (p in 1:N_p) {
    for (i in 1:N_i) {
      for (j in 1:N_j) {
        tau[p, i, j][1] = exp(log(tau_g[p][1]) + sig_tau_TRANS * tau_TRANS[p, i, j][1]);
        tau[p, i, j][2] = exp(log(tau_g[p][2]) + sig_tau_TRANS * tau_TRANS[p, i, j][2]);
      }
    }
  }
}

model {
  // Likelihood: observed data
  for (n in 1:N_obs) {
    int p = z_ind_p[n];
    int i = z_ind_i[n];
    int j = z_ind_j[n];

    z[n] ~ multi_normal(mu[p, i, j],
                        quad_form_diag(Omega[p, i, j], tau[p, i, j]));
  }

  // Likelihood: censored data 
  for (n in 1:N_cens) {
    int p = c_ind_p[n];
    int i = c_ind_i[n];
    int j = c_ind_j[n];

    target += normal_lcdf(alpha[p, i, j] | mu[p, i, j][1], tau[p, i, j][1]);
  }
  
  // Prior specification
  for (p in 1:N_p) {
    for (i in 1:N_i) {
      for (j in 1:N_j) {
        rho_TRANS[p, i, j] ~ normal(0, 1);
      }
    }
  }
  rho_g_TRANS ~ normal(0, .5);
  
  for (p in 1:N_p){
    for (i in 1:N_i){
      for(j in 1:N_j){
        mu[p, i, j] ~ multi_normal(m[p, j], diag_matrix(rep_vector(square(sig_mu), 2)));
        tau_TRANS[p, i, j] ~ normal(0, 1);
        alpha[p, i, j] ~ normal(alpha_g[i], sig_alpha);
      }
    }
  }
      
  for (p in 1:N_p) {
    tau_g[p][1] ~ lognormal(log(1.5), 0.75);   // x (depth): prior median 1.5m
    tau_g[p][2] ~ lognormal(log(1.0), 0.5);    // y (lateral): prior median 1.0m
  }
  for (i in 1:N_i) {
    alpha_g[i] ~ normal(4.0, 1.0);
  }
  
  // Hyperpriors
  sig_mu ~ exponential(5);
  sig_alpha ~ exponential(5);
  sig_rho_TRANS ~ exponential(5);
  sig_tau_TRANS ~ exponential(5);
}
