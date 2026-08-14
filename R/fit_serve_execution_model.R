# =============================================================================
# Fit the hierarchical serve-execution model (Section 4 of the paper).
#
# Reads  : data/serves.csv, data/targets.csv
# Writes : model_output/serve_execution_mcmc.rds  
# =============================================================================

library(cmdstanr)
library(tidyverse)

REGION_LEVELS <- c("deuce wide", "deuce tee", "ad tee", "ad wide")

serves  <- read_csv("data/serves.csv") |> 
  mutate(region = factor(region, levels = REGION_LEVELS))
targets <- read_csv("data/targets.csv") |> 
  mutate(region = factor(region, levels = REGION_LEVELS))

obs_dat  <- serves |> filter(clear_net == 1, let == 0)
cens_dat <- serves |> filter(clear_net == 0)

N_p <- length(unique(serves$player_id))
N_i <- 2
N_j <- length(unique(serves$region))

# Upper bound on the latent net-contact threshold alpha[p,i,j]: the shallowest
# serve actually seen to land in that cell.
a <- array(dim = c(N_p, N_i, N_j))
for (p in 1:N_p) {
  for (i in 1:N_i) {
    for (j in 1:N_j) {
      a[p, i, j] <- serves |>
        filter(player_id == p,
               serve_period == i,
               as.numeric(region) == j) |>
        slice_min(shot_x, with_ties = FALSE) |>
        pull(shot_x)
    }
  }
}

# Stated targets, used as the prior means on mu. No serve-period index: each
# player marked one target per region, used for both first and second serves.
m <- array(dim = c(N_p, N_j, 2))
for (p in 1:N_p) {
  for (j in 1:N_j) {
    m[p, j, ] <- targets |>
      filter(player_id == p, as.numeric(region) == j) |>
      select(target_x, target_y) |>
      unlist()
  }
}

stan_data <- list(N_obs = nrow(obs_dat),
                  N_cens = nrow(cens_dat),
                  N_p = N_p,
                  N_i = N_i,
                  N_j = N_j,
                  z = obs_dat |> select(shot_x, shot_y) |>
                    as.matrix(),
                  m = m,
                  a = a,
                  z_ind_p = obs_dat$player_id,
                  z_ind_i = obs_dat$serve_period,
                  z_ind_j = as.numeric(obs_dat$region),
                  c_ind_p = cens_dat$player_id,
                  c_ind_i = cens_dat$serve_period,
                  c_ind_j = as.numeric(cens_dat$region))

mod <- cmdstan_model("stan/serve_execution_model.stan")

fit <- mod$sample(
  data = stan_data,
  seed = 42,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2500,
  iter_sampling = 2500
)

fit$diagnostic_summary()
fit$summary()

fit$save_output_files(dir = "model_output/",
                      basename = "serve_execution_mcmc",
                      timestamp = FALSE)
fit$save_object("model_output/serve_execution_mcmc.rds")
