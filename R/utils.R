# =============================================================================
# Helper functions
# =============================================================================

# Court lines for the returner's half of the court, in metres, with the origin
# at the intersection of the centre service line and the net. Matches the
# coordinate frame of data/serves.csv and data/targets.csv.
geom_halfcourt <- function(my_col = "gray40") {
  court_dat <- data.frame(
    x    = c(0, 0, 11.887, 0, 0, 0, 0, 6.4),
    xend = c(11.887, 0, 11.887, 11.887, 11.887, 11.887, 6.4, 6.4),
    y    = c(5.486, 5.486, 5.486, -5.486, 4.115, -4.115, 0, 4.115),
    yend = c(5.486, -5.486, -5.486, -5.486, 4.115, -4.115, 0, -4.115)
  )
  ggplot2::geom_segment(
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
    data = court_dat, color = my_col
  )
}

# Execution-error parameters at an aiming location between the two anchors
# (Section 4.4). The Wide and T anchors are interpolated by the angle the aim
# subtends at the serve impact point; an aim outside the two anchors takes the
# nearer anchor's parameters unchanged.
interpolate_cov <- function(x, y,
                            mu_W, mu_T,
                            sig_W, sig_T,
                            corr_W, corr_T,
                            alpha_W, alpha_T,
                            impact_pos = c(-11.887, 0)) {

  angle <- function(mu) atan2(mu[2] - impact_pos[2], mu[1] - impact_pos[1])

  angle_aim <- atan2(y - impact_pos[2], x - impact_pos[1])
  angle_W   <- angle(mu_W)
  angle_T   <- angle(mu_T)

  # order the anchors so that angle_T < angle_W
  if (angle_T > angle_W) {
    tmp <- list(angle_T, sig_T, corr_T, alpha_T)
    angle_T <- angle_W;   sig_T <- sig_W;      corr_T <- corr_W;     alpha_T <- alpha_W
    angle_W <- tmp[[1]];  sig_W <- tmp[[2]];   corr_W <- tmp[[3]];   alpha_W <- tmp[[4]]
  }

  w <- (angle_aim - angle_T) / (angle_W - angle_T)
  w <- min(max(w, 0), 1)   # snap to the nearer anchor outside the span

  list(sig   = (1 - w) * sig_T   + w * sig_W,
       corr  = (1 - w) * corr_T  + w * corr_W,
       alpha = (1 - w) * alpha_T + w * alpha_W)
}


# The serve objective Q_i evaluated over a grid of candidate aiming locations,
# for one player, court, serve period and posterior draw (Section 6).
#
#   Q(a) = ss_val + SUM_cells f(z | a) (V(z) - ss_val) * cell_area,
#
# using V = 2 * wp - 1 and the fact that f integrates to one, so mass landing
# outside the service box picks up the fault value ss_val without being
# enumerated and `wp` need only cover the box. ss_val is -1 on the second
# serve, and the second serve's optimal value on the first.
#
#   aims           candidate aiming locations, columns x and y
#   wp             win-probability surface for this player, court and period,
#                  on a regular grid of spacing grid_res, columns x, y, wp
#   pars_W/pars_T  one row each of execution parameters at the Wide and T
#                  anchors: mu_x, mu_y, tau_x, tau_y, rho, alpha
#
# Returns `aims` with an added `ev` column.
solve_q_surface <- function(aims, wp, pars_W, pars_T, ss_val, grid_res = 0.1) {

  v          <- 2 * wp$wp - 1
  bounce_mat <- t(as.matrix(wp[, c("x", "y")]))
  cell_area  <- grid_res^2

  ev <- numeric(nrow(aims))
  for (k in seq_len(nrow(aims))) {
    ic <- interpolate_cov(aims$x[k], aims$y[k],
                          c(pars_W$mu_x,  pars_W$mu_y),
                          c(pars_T$mu_x,  pars_T$mu_y),
                          c(pars_W$tau_x, pars_W$tau_y),
                          c(pars_T$tau_x, pars_T$tau_y),
                          pars_W$rho,   pars_T$rho,
                          pars_W$alpha, pars_T$alpha)

    Sigma <- diag(ic$sig) %*%
             matrix(c(1, ic$corr, ic$corr, 1), 2) %*%
             diag(ic$sig)

    # A serve shallower than the net-clearance threshold alpha contacts the
    # net. alpha depends on the aim, so the mask is rebuilt for every aim, and
    # area-weighted across the grid column it straddles to keep Q smooth.
    in_play <- pmin(pmax((wp$x + grid_res / 2 - ic$alpha) / grid_res, 0), 1)
    dens    <- exp(mgcv::dmvn(bounce_mat,
                              mu = c(aims$x[k], aims$y[k]), V = Sigma))

    ev[k] <- ss_val + cell_area * sum(dens * in_play * (v - ss_val))
  }

  aims$ev <- ev
  aims
}


# Local maxima of an objective surface, highest first. A cell is a mode when it
# strictly exceeds all of its existing 8 neighbours, so a cell on the edge of
# the grid is compared against the neighbours it has. Near-ties are then
# collapsed greedily from the top down, keeping a mode only when it lies
# further than merge_radius from every mode already kept.
find_modes <- function(surface, grid_res = 0.1, merge_radius = 0.5) {

  v  <- surface$ev
  ix <- round(surface$x / grid_res)
  iy <- round(surface$y / grid_res)
  cell_value <- stats::setNames(v, paste(ix, iy))

  offsets <- expand.grid(dx = -1:1, dy = -1:1)
  offsets <- offsets[!(offsets$dx == 0 & offsets$dy == 0), ]

  is_mode <- rep(TRUE, length(v))
  for (k in seq_len(nrow(offsets))) {
    nb <- cell_value[paste(ix + offsets$dx[k], iy + offsets$dy[k])]
    is_mode <- is_mode & (is.na(nb) | nb < v)
  }

  modes <- surface[is_mode, , drop = FALSE]
  modes <- modes[order(-modes$ev), , drop = FALSE]

  if (merge_radius > 0 && nrow(modes) > 1) {
    keep <- 1L
    for (i in seq_len(nrow(modes))[-1]) {
      d <- sqrt((modes$x[keep] - modes$x[i])^2 + (modes$y[keep] - modes$y[i])^2)
      if (min(d) > merge_radius + 1e-9) keep <- c(keep, i)
    }
    modes <- modes[keep, , drop = FALSE]
  }

  rownames(modes) <- NULL
  modes
}
