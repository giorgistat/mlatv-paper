# Bivariate geostatistical latent variable models: analysis code

This repository reproduces the analysis pipeline of *Bivariate geostatistical
latent variable models for the analysis of antibody density data*
(Giorgi & Wallin) on a **simulated dataset** with the same structure as the
real RedHot serosurvey used in the paper. The real survey data are not
redistributable, so `R/00_simulate_data.R` generates a synthetic dataset
(`sim_redhot`) with the same columns, scale, and rough age/spatial structure,
so that every downstream script runs unmodified.

## Requirements

- R (>= 4.2)
- The `mlatv` package (TMB/C++ backend; not included here — install from its
  own repository, e.g. `remotes::install_github("giorgistat/mlatv")`)
- `fmesher`, `rSPDE`, `Matrix`, `ggplot2`

## Structure

```
R/
  00_simulate_data.R              simulate sim_redhot (replaces data("redhot"))
  util_example.r                  REAL shared helpers (location prep,
                                   geographic filter, mesh plot)
  predictive_hdr_spatial_mc_example.R
                                   REAL predictive-draw / HDR helpers
                                   (make_new_location_projection,
                                   draw_laplace_field_weights, etc.)
  utils_helpers.R                 repo-specific additions: resolve_repo_root(),
                                   load_sim_redhot(), and shared PPC plotting
                                   helpers; sources the two files above
  01_univariate_ama1.R            Section 5.1 single-antigen fit (AMA1)
  02_univariate_msp1.R            Section 5.2 single-antigen fit (MSP1)
  03_bivariate_nonspatial.R       Section 3 bivariate mixture, non-spatial
  04_bivariate_spatial_joint.R    Section 4/5.3 joint spatial-Z fit + Table
                                   `tab:bootstrap_redhot`-style parameter table
  05_spatial_field_maps.R         posterior mean/SD maps of S_1, S_2
                                   (Figure `fig:application-prediction`-style)
  06_cross_validation_hdr.R       Section 5.5 CRPS/SCRPS/ES/HDR cross-validation
                                   (Table `tab:cv_redhot`-style)
  07_simulation_study.R           Section 6 simulation study, adapted from the
                                   real sim_redhot_v2.R (see below)
run_all.R                          runs 00 through 06 in order (07 is separate;
                                   see its own header for why)
```

Each numbered script writes its fitted objects to `output/rds/` and any
figures to `output/plot/`, and can be re-run independently once the scripts
before it have been run at least once (each script checks for and loads the
previous script's `.rds` output rather than re-fitting from scratch).

## Provenance of the shared helper files

`R/util_example.r` and `R/predictive_hdr_spatial_mc_example.R` are the
**real, original files** (supplied after an earlier version of this repo,
which had to reconstruct their contents by inference — see git history /
prior notes if you want to compare). They are included verbatim and sourced
by `utils_helpers.R`; no call sites elsewhere in the repo needed to change,
since the reconstructions happened to match the real functions' contracts
closely enough for the pipeline to keep working once the real files dropped
in.

One consequence worth knowing: `example_geographic_filter()`'s default
thresholds (`drop_south_lat = -0.62`, `drop_nw_lon_max = 34.80`,
`drop_nw_lat_min = -0.46`) are calibrated to the real Kisumu HDSS study
geography, not to this repo's synthetic domain. `R/00_simulate_data.R`
deliberately anchors its simulated longitude/latitude east of that
exclusion box (`lon0 = 34.95`, not the real survey's `~34.75`) so the real
filter is a no-op on simulated data rather than silently dropping a large
fraction of it. If you change the simulator's domain, re-check this.

## `R/07_simulation_study.R`: adapted from the real `sim_redhot_v2.R`

The uploaded `sim_redhot_v2.R` is the real Section 6 simulation-study script.
Its loop logic, formulas, grid design, and output schema are ported here
essentially unchanged. The only substantive change is where the "real
domain" and "truth" fit come from:

| | Original | This repo |
|---|---|---|
| Domain (`spatial_data`) | `readRDS("article_data.rds")$spatial` | built from `04`'s saved `mesh`/`X`/`age`/`delta_age`/`loc_unique` |
| Reference ("truth") fit | `readRDS("spatial_joint_fit.rds")` | `04_bivariate_spatial_joint.rds`'s `fit_joint` |

The grid confirmed by the real script is `N_SAMPLE_GRID = c(300, 600, 1200)`,
`N_PER_LOC = 5`, `RHO_S_GRID`/`RHO_T_GRID = {0.1, 0.5, 0.9}`,
`SIGMA_SCALE_GRID = {1, 0.5, 0.25, 0.1}` — this uploaded version already has
the `N_SAMPLE_GRID` bug (`seq(300L, 600L, 1200L)`, which evaluates to just
`300`) fixed, consistent with prior notes that a re-run with the corrected
grid was still on the horizon. The real script also confirms the actual
`mlatv` function name and signature used for the non-spatial two-stage MLE
step — `fit_bivariate_lbm_latent_t_age_delta_mle_v3()`, with explicit
`breaks`/`mids` rather than a bare quadrature size — which `03` has been
updated to match (it previously guessed at
`fit_bivariate_lbm_latent_t_mle(..., M = ...)`).

The real script deliberately seeds nothing, so parallel/repeated runs (e.g.
SLURM array-job tasks) produce independent replicates without seed
coordination; that design is preserved. Because the full grid
(3×3×3×4 = 108 scenarios, each fitting three spatial models) is far too slow
for a local smoke test, `07` adds a `DEMO_MODE` toggle (on by default, not
part of the original) that shrinks every grid dimension to length 1. Set
`SIM07_FULL_GRID=true` in the environment to run the real published grid —
expect this to take a long time and to be best run as an HPC array job, one
task per (scenario, replicate), as in the original.



