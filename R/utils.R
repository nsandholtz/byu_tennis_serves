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
