# Purpose: Location-level five-fold spatial cross-validation comparing the
#   joint bivariate spatial model, two separate univariate spatial models,
#   and the non-spatial bivariate model, using CRPS, energy score, marginal
#   interval coverage/width, and the bivariate 95% highest-density region
#   (Section 5.5, Table `tab:cv_redhot`-style).
# Adapted from: example31_redhot_spatial_predictive_hdr_cv.R
# Changes from the original:
#   - Self-contained: does NOT source() article_code/predictive_hdr_spatial_mc_example.R
#     (not among the uploaded files -- see README "What's missing", item 3).
#   - Reads starting parameters from THIS repo's own 03/04 fits, not from
#     example19/example23 RDS files.
#   - Defaults to a small subsample (N_LOCATION) for a fast demo run.
# Output: output/rds/06_cv_hdr_summary.csv, output/rds/06_cv_hdr.rds

library(mlatv)
suppressPackageStartupMessages({ library(Matrix); library(fmesher) })
source(file.path(if (exists("resolve_repo_root")) resolve_repo_root() else getwd(),
                  "R", "utils_helpers.R"))
# utils_helpers.R sources both R/util_example.r and
# R/predictive_hdr_spatial_mc_example.R, which supply
# example_prepare_locations()/example_geographic_filter() and
# make_new_location_projection()/draw_laplace_field_weights() respectively.
# (Earlier versions of this repo had to guess at these; both files are now
# the real originals.)
repo_root <- resolve_repo_root()
rds_dir  <- file.path(repo_root, "output", "rds")

# ---- Run settings ---------------------------------------------------------
N_LOCATION   <- 150L
N_FOLD       <- 5L
N_DRAW       <- 500L
GRID_N       <- 40L
REGION_PROB  <- 0.95
M_GRID       <- 10L
MESH_CUTOFF  <- 100
MESH_MAX_EDGE <- c(250, 700)
SEED <- 20260824L
set.seed(SEED)

# ---- Scoring helpers (unchanged from example31.R) --------------------------
.crps_mc <- function(draws, observed) {
  draws <- as.numeric(draws)
  n <- length(draws)
  ordered <- sort(draws)
  coefficient <- 2 * seq_len(n) - n - 1
  mean(abs(draws - observed)) - sum(coefficient * ordered) / n^2
}
.energy_score_mc <- function(draws_1, draws_2, observed_1, observed_2) {
  draws <- cbind(draws_1, draws_2)
  first_term <- mean(sqrt((draws[, 1L] - observed_1)^2 + (draws[, 2L] - observed_2)^2))
  shifted <- draws[c(seq.int(2L, nrow(draws)), 1L), , drop = FALSE]
  second_term <- mean(sqrt(rowSums((draws - shifted)^2)))
  first_term - 0.5 * second_term
}
.prediction_summary <- function(y_1, y_2, draws_1, draws_2, model) {
  probability <- c(0.025, 0.975)
  interval_1 <- t(apply(draws_1, 1L, stats::quantile, probs = probability))
  interval_2 <- t(apply(draws_2, 1L, stats::quantile, probs = probability))
  data.frame(
    model = model, outcome = c("AMA1", "MSP1", "Pair"),
    crps = c(
      mean(vapply(seq_along(y_1), function(i) .crps_mc(draws_1[i, ], y_1[[i]]), numeric(1L))),
      mean(vapply(seq_along(y_2), function(i) .crps_mc(draws_2[i, ], y_2[[i]]), numeric(1L))),
      NA_real_
    ),
    energy_score = c(NA_real_, NA_real_,
      mean(vapply(seq_along(y_1), function(i)
        .energy_score_mc(draws_1[i, ], draws_2[i, ], y_1[[i]], y_2[[i]]), numeric(1L)))
    ),
    coverage = c(
      mean(y_1 >= interval_1[, 1L] & y_1 <= interval_1[, 2L]),
      mean(y_2 >= interval_2[, 1L] & y_2 <= interval_2[, 2L]),
      NA_real_
    ),
    width = c(
      mean(interval_1[, 2L] - interval_1[, 1L]),
      mean(interval_2[, 2L] - interval_2[, 1L]),
      NA_real_
    )
  )
}

# ---- Prepare data -----------------------------------------------------------
dat <- load_sim_redhot(repo_root)
dat <- dat[
  is.finite(dat$age) & dat$age > 1 &
    is.finite(dat$ama_norm) & dat$ama_norm > 0 &
    is.finite(dat$msp_norm) & dat$msp_norm > 0 &
    is.finite(dat$long_m) & is.finite(dat$lat_m),
  , drop = FALSE
]
dat$y_1 <- log(dat$ama_norm)
dat$y_2 <- log(dat$msp_norm)
site_center <- colMeans(dat[, c("long_m", "lat_m")])
dat$x_m <- dat$long_m - site_center[["long_m"]]
dat$y_m <- dat$lat_m  - site_center[["lat_m"]]
dat <- example_geographic_filter(dat)

prep <- example_prepare_locations(dat, n_sub = N_LOCATION, max_per_loc = 4L)
dat        <- prep$dat
loc        <- prep$loc
loc_unique <- prep$loc_unique
age        <- prep$age
y_1        <- prep$y_1
y_2        <- prep$y_2
n_obs      <- prep$n_obs
n_loc      <- prep$n_unique

age_knot <- 10
X <- cbind(
  `(Intercept)`   = 1,
  log_age         = log(age),
  log_age_above10 = pmax(log(age) - log(age_knot), 0)
)
delta_age <- log(age)

location_key <- paste(dat$x_m, dat$y_m, sep = "_")
location_id  <- match(location_key, unique(location_key))
fold_by_location <- sample(rep(seq_len(N_FOLD), length.out = n_loc))
fold <- fold_by_location[location_id]

mesh <- fmesher::fm_mesh_2d_inla(loc = loc_unique, max.edge = MESH_MAX_EDGE,
                                  cutoff = MESH_CUTOFF)
cat(sprintf("CV: n=%d, locations=%d, mesh nodes=%d\n", n_obs, n_loc, mesh$n))

# ---- Starting parameters, from THIS repo's own upstream fits ----------------
fit19_path <- file.path(rds_dir, "03_bivariate_nonspatial.rds")
fit23_path <- file.path(rds_dir, "04_bivariate_spatial_joint.rds")
if (!file.exists(fit19_path)) stop("Run 03_bivariate_nonspatial.R first.")
if (!file.exists(fit23_path)) stop("Run 04_bivariate_spatial_joint.R first.")
joint_nonspatial_ref <- readRDS(fit19_path)$fit
joint_spatial_ref    <- readRDS(fit23_path)$fit_fixed
joint_spatial_start <- as.list(joint_spatial_ref$par[
  c("range_1", "range_2", "rho", "beta_A", "beta_M")
])
univariate_spatial_start <- list(
  range = mean(joint_spatial_ref$par[c("range_1", "range_2")]),
  beta  = mean(joint_spatial_ref$par[c("beta_A", "beta_M")])
)

# ---- Fit the separate-univariate parameters once (fixed across folds) -------
cat("\nEstimating separate-univariate parameters once...\n")
full_separate_fit <- fit_bivariate_spatial_z_separate(
  y_1 = y_1, y_2 = y_2, mesh = mesh, loc = loc, X = X,
  nonspatial_method = "mle", nonspatial_convergence = "accept", M = M_GRID,
  nonspatial_control = list(eval.max = 150L, iter.max = 150L),
  spatial_starts = list(A = univariate_spatial_start, M = univariate_spatial_start),
  spatial_control = list(eval.max = 100L, iter.max = 100L), verbose = 0L
)

joint_draws_1 <- matrix(NA_real_, n_obs, N_DRAW)
joint_draws_2 <- matrix(NA_real_, n_obs, N_DRAW)
univariate_draws_1 <- matrix(NA_real_, n_obs, N_DRAW)
univariate_draws_2 <- matrix(NA_real_, n_obs, N_DRAW)
nonspatial_draws_1 <- matrix(NA_real_, n_obs, N_DRAW)
nonspatial_draws_2 <- matrix(NA_real_, n_obs, N_DRAW)

for (fold_index in seq_len(N_FOLD)) {
  cat(sprintf("\n--- Fold %d/%d ---\n", fold_index, N_FOLD))
  test_index  <- which(fold == fold_index)
  train_index <- which(fold != fold_index)

  joint_nonspatial <- joint_nonspatial_ref
  joint_nonspatial$X <- X[train_index, , drop = FALSE]
  joint_nonspatial$delta_age <- delta_age[train_index]

  joint_fit <- fit_bivariate_spatial_z_joint(
    y_1 = y_1[train_index], y_2 = y_2[train_index],
    X = X[train_index, , drop = FALSE], delta_age = delta_age[train_index],
    mesh = mesh, loc = loc[train_index, , drop = FALSE], M = M_GRID,
    start = joint_spatial_start, estimate_theta0 = FALSE,
    nonspatial_fit = joint_nonspatial, nonspatial_convergence = "accept",
    control = list(eval.max = 1L, iter.max = 0L), verbose = 0L
  )

  separate_fit <- fit_bivariate_spatial_z_separate(
    y_1 = y_1[train_index], y_2 = y_2[train_index], mesh = mesh,
    loc = loc[train_index, , drop = FALSE], X = X[train_index, , drop = FALSE],
    nonspatial_fits = list(A = full_separate_fit$fit_A$nonspatial_fit,
                            M = full_separate_fit$fit_M$nonspatial_fit),
    nonspatial_method = "mle", nonspatial_convergence = "accept", M = M_GRID,
    spatial_starts = list(
      A = c(as.list(full_separate_fit$fit_A$par[c("range", "beta")]),
            list(w = full_separate_fit$fit_A$w_hat)),
      M = c(as.list(full_separate_fit$fit_M$par[c("range", "beta")]),
            list(w = full_separate_fit$fit_M$w_hat))
    ),
    spatial_control = list(eval.max = 1L, iter.max = 0L), verbose = 0L
  )

  # Predictive draws at the held-out locations
  A_joint <- make_new_location_projection(joint_fit$bsz_data, loc[test_index, , drop = FALSE])
  joint_weight <- draw_laplace_field_weights(joint_fit, n_draws = N_DRAW, field = "laplace")
  n_mesh <- joint_fit$bsz_data$n_mesh
  field_1 <- as.matrix(A_joint %*% joint_weight[seq_len(n_mesh), , drop = FALSE])
  field_2 <- as.matrix(A_joint %*% joint_weight[n_mesh + seq_len(n_mesh), , drop = FALSE])

  A_A <- make_new_location_projection(separate_fit$fit_A$usz_data, loc[test_index, , drop = FALSE])
  A_M <- make_new_location_projection(separate_fit$fit_M$usz_data, loc[test_index, , drop = FALSE])
  field_A <- as.matrix(A_A %*% draw_laplace_field_weights(
    separate_fit$fit_A, n_draws = N_DRAW, field = "laplace"))
  field_M <- as.matrix(A_M %*% draw_laplace_field_weights(
    separate_fit$fit_M, n_draws = N_DRAW, field = "laplace"))

  for (k in seq_along(test_index)) {
    i <- test_index[[k]]
    X_mc <- X[rep(i, N_DRAW), , drop = FALSE]
    sim_joint <- mlatv:::.simulate_bivariate_lbm_latent_t(
      par = joint_nonspatial$par, X = X_mc,
      delta_age = rep(delta_age[[i]], N_DRAW),
      eta_extra_1 = joint_fit$par[["beta_A"]] * field_1[k, ],
      eta_extra_2 = joint_fit$par[["beta_M"]] * field_2[k, ],
      breaks = joint_nonspatial$breaks, mids = joint_nonspatial$mids,
      latent_method = "quadrature"
    )
    joint_draws_1[i, ] <- sim_joint$y_1
    joint_draws_2[i, ] <- sim_joint$y_2

    univariate_draws_1[i, ] <- mlatv:::.simulate_univariate_lbm_latent_t(
      par = separate_fit$fit_A$nonspatial_fit$par, X = X_mc,
      eta_extra = separate_fit$fit_A$par[["beta"]] * field_A[k, ]
    )$y
    univariate_draws_2[i, ] <- mlatv:::.simulate_univariate_lbm_latent_t(
      par = separate_fit$fit_M$nonspatial_fit$par, X = X_mc,
      eta_extra = separate_fit$fit_M$par[["beta"]] * field_M[k, ]
    )$y
  }

  X_ns <- X[rep(test_index, each = N_DRAW), , drop = FALSE]
  sim_ns <- mlatv:::.simulate_bivariate_lbm_latent_t(
    par = joint_nonspatial$par, X = X_ns,
    delta_age = rep(delta_age[test_index], each = N_DRAW),
    breaks = joint_nonspatial$breaks, mids = joint_nonspatial$mids,
    latent_method = "quadrature"
  )
  nonspatial_draws_1[test_index, ] <- matrix(sim_ns$y_1, nrow = length(test_index),
                                              ncol = N_DRAW, byrow = TRUE)
  nonspatial_draws_2[test_index, ] <- matrix(sim_ns$y_2, nrow = length(test_index),
                                              ncol = N_DRAW, byrow = TRUE)
  cat(sprintf("Fold %d done.\n", fold_index))
}

# ---- Score and HDR summaries -------------------------------------------------
predictive_summary <- rbind(
  .prediction_summary(y_1, y_2, joint_draws_1, joint_draws_2, "Joint spatial"),
  .prediction_summary(y_1, y_2, univariate_draws_1, univariate_draws_2, "Univariate spatial"),
  .prediction_summary(y_1, y_2, nonspatial_draws_1, nonspatial_draws_2, "Non-spatial")
)

joint_hdr       <- evaluate_bivariate_predictive_regions(
  y_1 = y_1, y_2 = y_2, draws_1 = joint_draws_1, draws_2 = joint_draws_2,
  prob = REGION_PROB, grid_n = GRID_N)
univariate_hdr  <- evaluate_bivariate_predictive_regions(
  y_1 = y_1, y_2 = y_2, draws_1 = univariate_draws_1, draws_2 = univariate_draws_2,
  prob = REGION_PROB, grid_n = GRID_N)
nonspatial_hdr  <- evaluate_bivariate_predictive_regions(
  y_1 = y_1, y_2 = y_2, draws_1 = nonspatial_draws_1, draws_2 = nonspatial_draws_2,
  prob = REGION_PROB, grid_n = GRID_N)

hdr_summary <- rbind(
  data.frame(model = "Joint spatial", joint_hdr$summary),
  data.frame(model = "Univariate spatial", univariate_hdr$summary),
  data.frame(model = "Non-spatial", nonspatial_hdr$summary)
)
hdr_summary$coverage_percent <- 100 * hdr_summary$coverage

cat("\nPredictive score and marginal interval summary\n")
print(predictive_summary, row.names = FALSE, digits = 4)
cat("\nJoint HDR summary\n")
print(hdr_summary, row.names = FALSE, digits = 4)

# ---- Save ---------------------------------------------------------------------
utils::write.csv(predictive_summary, file.path(rds_dir, "06_cv_scores_summary.csv"),
                  row.names = FALSE)
utils::write.csv(hdr_summary, file.path(rds_dir, "06_cv_hdr_summary.csv"),
                  row.names = FALSE)
saveRDS(
  list(predictive_summary = predictive_summary, hdr_summary = hdr_summary,
       fold = fold, settings = list(N_LOCATION = N_LOCATION, N_FOLD = N_FOLD,
                                     N_DRAW = N_DRAW, GRID_N = GRID_N,
                                     REGION_PROB = REGION_PROB, seed = SEED)),
  file.path(rds_dir, "06_cv_hdr.rds")
)
cat("\nSaved: output/rds/06_cv_scores_summary.csv\n")
cat("Saved: output/rds/06_cv_hdr_summary.csv\n")
cat("Saved: output/rds/06_cv_hdr.rds\n")
