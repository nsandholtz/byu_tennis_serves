# =============================================================================
# Illustration of the fitted serve-execution model, Player 2 (Figure 4).
#
# Run from the repository root:
#   Rscript R/illustrate_exec_err.R
#
# For each serve period and target region, draws the 90% contour of the fitted
# bivariate normal execution distribution at the posterior mean, together with
# the observed bounce locations and the stated target. The portion of each
# contour lying short of the fitted censoring threshold alpha[p,i,j] is hatched:
# those serves contact the net and produce no observable bounce.
#
# Reads  : model_output/posterior_draws_1000.csv (written by
#          R/fit_serve_execution_model.R), data/serves.csv, data/targets.csv
# =============================================================================

source("R/utils.R")   

suppressMessages({
  library(tidyverse)
  library(ggpattern)
  library(ggtext)
})

PLAYER <- 2
REGION_LEVELS <- c("deuce wide", "deuce tee", "ad tee", "ad wide")

draws   <- read_csv("model_output/posterior_draws_1000.csv", show_col_types = FALSE)
serves  <- read_csv("data/serves.csv",  show_col_types = FALSE) |>
  mutate(region = factor(region, levels = REGION_LEVELS))
targets <- read_csv("data/targets.csv", show_col_types = FALSE) |>
  mutate(region = factor(region, levels = REGION_LEVELS))

# ---- posterior means of the cell-level parameters ---------------------------

pars <- draws |>
  filter(player_id == PLAYER) |>
  group_by(player_id, serve_period, region_i = region) |>
  summarise(across(c(mu_x, mu_y, tau_x, tau_y, rho, alpha), mean),
            .groups = "drop")

# ---- 90% contour of each fitted execution distribution ----------------------
# (z - mu)' Sigma^{-1} (z - mu) = qchisq(.90, 2), traced by mapping the unit
# circle through the Cholesky factor. Avoids a dependency on the ellipse package.
ellipse_pts <- function(mu_x, mu_y, tau_x, tau_y, rho, level = 0.90, n = 240) {
  Sigma <- matrix(c(tau_x^2, rho * tau_x * tau_y,
                    rho * tau_x * tau_y, tau_y^2), nrow = 2)
  th <- seq(0, 2 * pi, length.out = n)
  pts <- t(c(mu_x, mu_y) + sqrt(qchisq(level, 2)) *
             (t(chol(Sigma)) %*% rbind(cos(th), sin(th))))
  tibble(x = pts[, 1], y = pts[, 2])
}

ellipse_df <- pars |>
  mutate(cell = row_number()) |>
  rowwise() |>
  reframe(cell = cell, serve_period = serve_period, region_i = region_i,
          alpha = alpha,
          ellipse_pts(mu_x, mu_y, tau_x, tau_y, rho)) |>
  mutate(region = factor(REGION_LEVELS[region_i], levels = REGION_LEVELS),
         first_serve_nice = if_else(serve_period == 1, "1st\nServe", "2nd\nServe"),
         server_name = paste("Player", PLAYER))

# The contour is split at the censoring threshold: mass short of alpha in the
# depth (x) direction hits the net.
ellipse_obs  <- ellipse_df |> filter(x >  alpha)
ellipse_cens <- ellipse_df |> filter(x <= alpha)

# ---- observed serves, formatted as in Figure 3 ------------------------------
serves <- serves |> filter(player_id == PLAYER, let == 0)

serves$first_serve_nice <- if_else(serves$serve_period == 1, "1st\nServe", "2nd\nServe")
serves$server_name <- paste("Player", serves$player_id)

# ---------------------------------------------------------------------------
# NOTE: the coordinates assigned below are INVENTED FOR DISPLAY ONLY.
#
# Censored serves (clear_net == 0) contacted the net, so no bounce was ever
# recorded: shot_x and shot_y are NA for all of them. For visual purposes, we
# display these net shots just short of the net. The anchors are chosen to sit
# near the lateral position where that region's successful serves land, so the
# color coding stays readable. Nothing about these positions is a measurement.
# Do not read anything into where a censored point falls beyond which region
# it was aimed at.

serves$x_serve_bounce <- serves$shot_x
serves$x_serve_bounce[serves$clear_net == 0] <- runif(sum(serves$clear_net == 0), -.2, 0)

serves$y_serve_bounce <- serves$shot_y
for (rg in REGION_LEVELS) {
  yy <- c("deuce wide" = 4.115 * (2 / 3), "deuce tee" = 4.115 * (1 / 10),
          "ad tee" = -4.115 * (1 / 10), "ad wide" = -4.115 * (2 / 3))[[rg]]
  k <- serves$clear_net == 0 & serves$region == rg
  serves$y_serve_bounce[k] <- yy + rnorm(sum(k), 0, .5)
}

targets <- targets |> filter(player_id == PLAYER) |>
  mutate(server_name = paste("Player", player_id))

PT_LEVELS <- c("Observed", "Censored", "Stated target")
serves$point_type  <- factor(if_else(serves$clear_net == 1, "Observed", "Censored"),
                             levels = PT_LEVELS)
targets$point_type <- factor("Stated target", levels = PT_LEVELS)

# ---- plot -------------------------------------------------------------------

p <- ggplot(mapping = aes(x = x, y = y)) +
  geom_halfcourt() +
  # fitted execution distribution: solid where the serve lands, hatched where
  # it contacts the net
  geom_polygon(data = ellipse_obs,
               aes(group = cell, fill = region), alpha = 0.3) +
  geom_polygon_pattern(data = ellipse_cens,
                       aes(group = cell, fill = region),
                       alpha = 0.3,
                       pattern = "stripe", pattern_angle = 45,
                       pattern_density = 0.35, pattern_spacing = 0.015,
                       pattern_fill = "white", pattern_colour = NA) +
  geom_path(data = ellipse_df, aes(group = cell, colour = region),
            linewidth = 0.18, alpha = 0.55) +
  geom_halfcourt() +
  # observed data, matching Figure 3
  geom_point(data = serves,
             aes(x = x_serve_bounce, y = y_serve_bounce,
                 colour = region, shape = point_type),
             alpha = .5, cex = 1) +
  geom_point(data = targets,
             aes(x = target_x, y = target_y, colour = region, shape = point_type),
             size = 1, stroke = 1) +
  facet_grid(server_name ~ first_serve_nice, switch = "y") +
  coord_equal(ylim = c(-5.6, 5.7)) +
  xlab("x (meters)") + ylab("y (meters)") +
  scale_x_continuous(breaks = seq(0, 12, 3)) +
  scale_y_continuous(breaks = seq(-6, 6, 3)) +
  scale_color_manual(values = palette.colors(),
                     name = "Target Region",
                     labels = c("deuce wide" = "Deuce Wide",
                                "deuce tee"  = "Deuce T",
                                "ad tee"     = "Ad T",
                                "ad wide"    = "Ad Wide")) +
  scale_fill_manual(values = palette.colors(), guide = "none") +
  scale_shape_manual(values = c("Observed" = 16, "Censored" = 1,
                                "Stated target" = 4)) +
  guides(
    color = guide_legend(title.position = "top",
                         override.aes = list(size = 2, alpha = 1, shape = 16)),
    shape = guide_legend(title.position = "top",
                         override.aes = list(size = 2, alpha = 1,
                                             shape = c(16, 1, 4)))
  ) +
  labs(shape = "Symbol") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        legend.position = "right",
        strip.text.y.left = element_blank(),
        legend.title = element_text(size = 6, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 6),
        strip.text = ggtext::element_markdown(size = 8, lineheight = 1.1),
        legend.background = element_rect(fill = "gray95", color = NA))

p
