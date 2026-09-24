# =============================================================================
# Solve the serve-aiming problem and locate the modes of each objective surface
# (Section 6 of the paper).
#
# Run from the repository root:
#   Rscript R/solve_modes.R
#
# For each player, posterior draw and service court, the objective is evaluated
# over a grid of candidate aiming locations by backward induction: the second
# serve first, where a fault loses the point, then the first serve, where a
# fault carries the value of the second-serve optimum. The modes of each
# surface are the candidate optimal aiming locations.
#
# Costs about two seconds per player and posterior draw, so roughly an hour for
# all eight players with the player loop spread across cores.
#
# Reads  : model_output/posterior_draws_1000.csv (written by
#          R/fit_serve_execution_model.R)
#          model_output/wp_surface_deuce.csv, model_output/wp_surface_ad.csv
# Writes : model_output/modes_1000.csv
# =============================================================================

source("R/utils.R")   # interpolate_cov(), solve_q_surface(), find_modes()

suppressMessages({
  library(tidyverse)
  library(parallel)
})

MERGE_RADIUS <- 0.5   # modes closer than this to a higher one are dropped
GRID_RES     <- 0.1   # spacing of both the aim grid and the win-probability grid

draws <- read_csv("model_output/posterior_draws_1000.csv", show_col_types = FALSE)

wp <- bind_rows(
  read_csv("model_output/wp_surface_deuce.csv", show_col_types = FALSE),
  read_csv("model_output/wp_surface_ad.csv",    show_col_types = FALSE)
)

# ---- aim grid ---------------------------------------------------------------
# Candidate aiming locations, generous around the service box and well outside
# it on every side, so that a mode can fall outside the box without the edge of
# the grid binding. The objective is flat at the fault value further out, so
# searching the whole half court would only add work.

aim_grid <- function(court) {
  y_range <- if (court == "Deuce") c(-1.7, 5.7) else c(-5.7, 1.7)
  expand_grid(x = seq(1.6, 8.2, by = GRID_RES),
              y = seq(y_range[1], y_range[2], by = GRID_RES))
}

AIMS <- list(Deuce = aim_grid("Deuce"), Ad = aim_grid("Ad"))

# ---- modes, every player, draw and court ------------------------------------

solve_player <- function(p) {
  draws_p <- filter(draws, player_id == p)
  out <- list()

  for (d in unique(draws_p$draw)) {
    pars <- filter(draws_p, draw == d)

    for (crt in c("Deuce", "Ad")) {
      anchor <- function(period, strat) {
        filter(pars, serve_period == period, court == crt, strategy == strat)
      }
      surface <- filter(wp, player_id == p, court == crt)

      # second serve: a fault ends the point, so it is worth -1
      q2 <- solve_q_surface(AIMS[[crt]],
                            filter(surface, serve_period == 2),
                            anchor(2, "Wide"), anchor(2, "T"),
                            ss_val = -1, grid_res = GRID_RES)

      # first serve: a fault carries the value of the second-serve optimum
      q1 <- solve_q_surface(AIMS[[crt]],
                            filter(surface, serve_period == 1),
                            anchor(1, "Wide"), anchor(1, "T"),
                            ss_val = max(q2$ev), grid_res = GRID_RES)

      out[[length(out) + 1]] <- bind_rows(
        find_modes(q1, GRID_RES, MERGE_RADIUS) |> mutate(serve_period = 1L),
        find_modes(q2, GRID_RES, MERGE_RADIUS) |> mutate(serve_period = 2L)
      ) |>
        mutate(player_id = p, court = crt, draw = d)
    }
  }

  message("player ", p, " done")
  bind_rows(out)
}

players <- sort(unique(draws$player_id))
res <- mclapply(players, solve_player,
                mc.cores = min(length(players), detectCores() - 1))

failed <- !map_lgl(res, is.data.frame)
if (any(failed)) stop("solve failed for players: ",
                      str_c(players[failed], collapse = ", "))

modes <- bind_rows(res) |>
  select(player_id, draw, court, serve_period, x, y, ev)

write_csv(modes, "model_output/modes_1000.csv")
