# Monte Carlo predictive HDRs for held-out spatial observations
#
# The fit objects supplied here must be fitted without the held-out location.
# Setting field = "laplace" integrates over the fitted latent field under the
# Gaussian Laplace approximation. Setting field = "mode" conditions on w_hat.

make_new_location_projection <- function(spatial_data, loc_new) {
  loc_new <- as.matrix(loc_new)
  if (!is.numeric(loc_new) || nrow(loc_new) < 1L ||
      ncol(loc_new) != spatial_data$dimension ||
      any(!is.finite(loc_new))) {
    stop(
      paste0(
        "loc_new must contain one or more finite locations with ",
        "spatial_data$dimension columns."
      ),
      call. = FALSE
    )
  }

  operator <- mlatv:::.build_rspde_operator(
    mesh = spatial_data$mesh,
    range = mlatv:::.default_mesh_range(spatial_data$mesh),
    nu = spatial_data$nu
  )
  methods::as(
    rSPDE::make_A(operator, loc = loc_new),
    "CsparseMatrix"
  )
}

draw_laplace_field_weights <- function(
    spatial_fit, n_draws, field = c("laplace", "mode")
) {
  field <- match.arg(field)
  if (!isTRUE(spatial_fit$laplace_valid) ||
      is.null(spatial_fit$w_hat)) {
    stop("spatial_fit must contain a valid Laplace mode.", call. = FALSE)
  }

  if (field == "mode") {
    return(matrix(
      rep(spatial_fit$w_hat, n_draws),
      nrow = length(spatial_fit$w_hat),
      ncol = n_draws
    ))
  }
  if (is.null(spatial_fit$H_hat)) {
    stop(
      "spatial_fit$H_hat is required for Laplace field draws.",
      call. = FALSE
    )
  }

  deviations <- sample_gmrf(spatial_fit$H_hat, nsim = n_draws)
  sweep(deviations, 1L, spatial_fit$w_hat, FUN = "+")
}

draw_bivariate_nonspatial_predictive <- function(
    nonspatial_fit, X_new, delta_age_new = NULL,
    n_draws = 4000L, seed = 20260816L
) {
  X_new <- as.matrix(X_new)
  if (nrow(X_new) != 1L) {
    stop("X_new must contain one held-out design row.", call. = FALSE)
  }
  X_mc <- X_new[rep(1L, n_draws), , drop = FALSE]
  delta_mc <- if (is.null(delta_age_new)) {
    NULL
  } else {
    rep(as.numeric(delta_age_new), n_draws)
  }

  set.seed(seed)
  simulation <- mlatv:::.simulate_bivariate_lbm_latent_t(
    par = nonspatial_fit$par,
    X = X_mc,
    delta_age = delta_mc,
    breaks = nonspatial_fit$breaks,
    mids = nonspatial_fit$mids,
    latent_method = "quadrature"
  )
  cbind(y_1 = simulation$y_1, y_2 = simulation$y_2)
}

draw_bivariate_spatial_predictive <- function(
    nonspatial_fit, spatial_fit, X_new, delta_age_new, loc_new,
    n_draws = 4000L, field = c("laplace", "mode"),
    seed = 20260816L
) {
  field <- match.arg(field)
  X_new <- as.matrix(X_new)
  if (nrow(X_new) != 1L) {
    stop("X_new must contain one held-out design row.", call. = FALSE)
  }

  set.seed(seed)
  A_new <- make_new_location_projection(spatial_fit$bsz_data, loc_new)
  w_mc <- draw_laplace_field_weights(
    spatial_fit, n_draws = n_draws, field = field
  )
  n_mesh <- spatial_fit$bsz_data$n_mesh

  if (identical(spatial_fit$field_structure, "shared")) {
    field_1 <- field_2 <- as.numeric(A_new %*% w_mc)
  } else {
    first <- seq_len(n_mesh)
    second <- n_mesh + first
    field_1 <- as.numeric(A_new %*% w_mc[first, , drop = FALSE])
    field_2 <- as.numeric(A_new %*% w_mc[second, , drop = FALSE])
  }

  X_mc <- X_new[rep(1L, n_draws), , drop = FALSE]
  simulation <- mlatv:::.simulate_bivariate_lbm_latent_t(
    par = nonspatial_fit$par,
    X = X_mc,
    delta_age = rep(as.numeric(delta_age_new), n_draws),
    eta_extra_1 = spatial_fit$par[["beta_A"]] * field_1,
    eta_extra_2 = spatial_fit$par[["beta_M"]] * field_2,
    breaks = nonspatial_fit$breaks,
    mids = nonspatial_fit$mids,
    latent_method = "quadrature"
  )
  cbind(y_1 = simulation$y_1, y_2 = simulation$y_2)
}

draw_univariate_spatial_predictive <- function(
    nonspatial_fit, spatial_fit, X_new, loc_new,
    n_draws = 4000L, field = c("laplace", "mode"),
    seed = 20260816L
) {
  field <- match.arg(field)
  X_new <- as.matrix(X_new)
  if (nrow(X_new) != 1L) {
    stop("X_new must contain one held-out design row.", call. = FALSE)
  }

  set.seed(seed)
  A_new <- make_new_location_projection(spatial_fit$usz_data, loc_new)
  w_mc <- draw_laplace_field_weights(
    spatial_fit, n_draws = n_draws, field = field
  )
  field_mc <- as.numeric(A_new %*% w_mc)
  X_mc <- X_new[rep(1L, n_draws), , drop = FALSE]
  simulation <- mlatv:::.simulate_univariate_lbm_latent_t(
    par = nonspatial_fit$par,
    X = X_mc,
    eta_extra = spatial_fit$par[["beta"]] * field_mc
  )
  simulation$y
}

# Example for one held-out individual -----------------------------------------
#
# The following object names are placeholders for fold-specific fits and the
# held-out person's covariates, location, and response. The entire location,
# not only this person, must have been omitted when fitting these objects.
if (FALSE) {
  joint_draws <- draw_bivariate_spatial_predictive(
    nonspatial_fit = fold_joint_nonspatial_fit,
    spatial_fit = fold_joint_spatial_fit,
    X_new = X_test[i, , drop = FALSE],
    delta_age_new = delta_age_test[[i]],
    loc_new = loc_test[i, , drop = FALSE],
    n_draws = 4000L,
    field = "laplace",
    seed = 1000L + i
  )
  joint_region <- estimate_bivariate_predictive_region(
    draws = joint_draws,
    observed = c(y_1_test[[i]], y_2_test[[i]]),
    prob = 0.95,
    grid_n = 120L
  )

  independent_draws <- cbind(
    y_1 = draw_univariate_spatial_predictive(
      nonspatial_fit = fold_separate_fit$fit_A$nonspatial_fit,
      spatial_fit = fold_separate_fit$fit_A,
      X_new = X_test[i, , drop = FALSE],
      loc_new = loc_test[i, , drop = FALSE],
      n_draws = 4000L,
      field = "laplace",
      seed = 2000L + i
    ),
    y_2 = draw_univariate_spatial_predictive(
      nonspatial_fit = fold_separate_fit$fit_M$nonspatial_fit,
      spatial_fit = fold_separate_fit$fit_M,
      X_new = X_test[i, , drop = FALSE],
      loc_new = loc_test[i, , drop = FALSE],
      n_draws = 4000L,
      field = "laplace",
      seed = 3000L + i
    )
  )
  independent_region <- estimate_bivariate_predictive_region(
    draws = independent_draws,
    observed = c(y_1_test[[i]], y_2_test[[i]]),
    prob = 0.95,
    grid_n = 120L
  )

  rbind(
    joint = c(
      covered = joint_region$covered,
      area = joint_region$area
    ),
    independent = c(
      covered = independent_region$covered,
      area = independent_region$area
    )
  )
}

# Cross-validation aggregation ------------------------------------------------
#
# Repeat the one-person calculation for every held-out individual. Store the
# joint draws in two n_test by n_draws matrices, retaining the pairing of the
# columns. Do the same for the two independent univariate predictive draws.
if (FALSE) {
  joint_evaluation <- evaluate_bivariate_predictive_regions(
    y_1 = y_1_test,
    y_2 = y_2_test,
    draws_1 = joint_draws_1,
    draws_2 = joint_draws_2,
    prob = 0.95,
    grid_n = 120L
  )
  independent_evaluation <- evaluate_bivariate_predictive_regions(
    y_1 = y_1_test,
    y_2 = y_2_test,
    draws_1 = independent_draws_1,
    draws_2 = independent_draws_2,
    prob = 0.95,
    grid_n = 120L
  )

  rbind(
    joint = joint_evaluation$summary,
    independent = independent_evaluation$summary
  )
}
