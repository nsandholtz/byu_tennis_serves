# byu_tennis_serves

Data and code to reproduce the results in "Stated, Realized, and Optimal Aiming
Strategies for the Tennis Serve: Experimental Evidence from Collegiate
Athletes," by Nathan Sandholtz, Ryan Hanson, Ron Hager, Stephanie Kovalchik,
and Gilbert Fellingham.

Eight members of the BYU men's tennis team each marked a target location in
four regions of the service court, then served at those targets. The paper
compares three aiming locations for every player, serve period and target
region: the target the player stated, the center of the serve distribution they
actually realized, and the target that maximizes their probability of winning
the point under their own execution error.

#### Running the code

All scripts are run from the repository root:

```
Rscript R/fit_serve_execution_model.R   # fit the execution model
Rscript R/solve_modes.R                 # optimal aiming locations, Section 6
```

`fit_serve_execution_model.R` comes first: it writes
`model_output/posterior_draws_1000.csv`, which `solve_modes.R` reads. That file is
included here, so the mode solving script can be run without refitting.
`solve_modes.R` takes about an hour, spread across cores by player; the other
two run in seconds.

`R/utils.R` holds the shared helpers: the court geometry used in the figures,
and the three functions that implement the aiming problem
(`interpolate_cov`, `solve_q_surface`, `find_modes`).

Requires R with `tidyverse`, `mgcv`, `ggpattern` and `ggtext`, and, for
refitting only, `cmdstanr` with a working CmdStan installation.
`stan/serve_execution_model` is a compiled binary and will need to be rebuilt
on another platform; `cmdstan_model()` does this automatically.

#### Data

`data/serves.csv`, one row per serve attempt. 567 serves from 8 players,
collected across eight sessions.

| column | description |
| --- | --- |
| `session_id` | collection session, 1 to 8 in chronological order |
| `player_id` | player, 1 to 8 |
| `serve_period` | 1 for first serve, 2 for second serve |
| `region` | target region: `deuce wide`, `deuce tee`, `ad tee`, `ad wide` |
| `index` | serve attempt number within the session |
| `speed_mph` | serve speed, miles per hour |
| `clear_net` | 1 if the serve cleared the net, 0 if it contacted the net |
| `let` | 1 if the serve was a let |
| `shot_x`, `shot_y` | bounce location, meters |

`data/targets.csv`, one row per player and region, giving the target the player
marked on the court (`target_x`, `target_y`) and their dominant hand. Each
player marked one target per region, used for both serve periods. 

Coordinates are in meters in the returner's half of the court, with the origin
at the intersection of the center service line and the net: `x` is depth from
the net and `y` is lateral position, positive toward the Deuce court. The same
frame is used throughout the model output. `geom_halfcourt()` in `R/utils.R`
draws the court lines in this frame.

Two notes on `serves.csv`:

- The 111 serves with `clear_net == 0` hit the net and did not cross over, so no bounce was ever
  recorded and `shot_x` and `shot_y` are empty. 
- The 19 serves with `let == 1` are excluded when fitting.

Region order is `deuce wide`, `deuce tee`, `ad tee`, `ad wide` throughout. Where
the model output codes region as an integer, 1 to 4 follows that order.

#### Model output

| file | description |
| --- | --- |
| `serve_execution_mcmc-{1..4}-390456.csv` | raw CmdStan output, four chains |
| `posterior_draws_1000.csv` | 1,000 posterior draws of the cell-level execution parameters (`mu_x`, `mu_y`, `tau_x`, `tau_y`, `rho`, `alpha`), one row per draw, player, serve period and region |
| `wp_surface_deuce.csv`, `wp_surface_ad.csv` | probability the server wins the point given a serve landing at each location, by player and serve period |
| `modes_1000.csv` | modes of the objective surface for each player, posterior draw, court and serve period, written by `R/solve_modes.R` |

The win-probability surfaces are derived from the shot-level trajectory and
win-probability model of Kovalchik, Ingram, Weeratunga and Goncu (2020),
"Space-Time VON CRAMM"
([arXiv:2005.12853](https://arxiv.org/abs/2005.12853)), applied recursively
over simulated rallies to give a point-level win probability for a serve
landing at each location on a 0.1 m grid over the service box, conditioned on
each player's average first- and second-serve speed. Section 6.1 and the
appendix of the paper describe the procedure.

#### License

MIT, see `LICENSE`.
