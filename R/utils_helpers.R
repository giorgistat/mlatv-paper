# Purpose: Shared helper functions used by scripts 01-06, plus repo-specific
#   additions not present in the original files.
#
# util_example.r and predictive_hdr_spatial_mc_example.R are now the REAL
# files (supplied after the initial version of this repo, which had to
# reconstruct example_prepare_locations(), example_geographic_filter(), and
# example_save_mesh_plot() from how they were called elsewhere). They are
# sourced below and used as-is. The functions defined directly in this file
# are repo-specific additions with no equivalent in the original scripts:
# resolve_repo_root(), load_sim_redhot(), and the posterior-predictive-check
# plotting helpers factored out of example28.R/example29.R.

.utils_helpers_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE)),
  error = function(e) file.path(getwd(), "R")
)
source(file.path(.utils_helpers_dir, "util_example.r"))
source(file.path(.utils_helpers_dir, "predictive_hdr_spatial_mc_example.R"))
#' `Rscript path/to/script.R` or `source()`d from an interactive session at
#' the repo root.
resolve_repo_root <- function() {
  repo_root <- getwd()
  if (!file.exists(file.path(repo_root, "R", "00_simulate_data.R"))) {
    script_path <- tryCatch(
      normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
      error = function(e) NA_character_
    )
    if (!is.na(script_path)) {
      repo_root <- normalizePath(
        file.path(dirname(script_path), ".."),
        winslash = "/", mustWork = TRUE
      )
    }
  }
  repo_root
}

#' Load the simulated dataset, generating it first if needed.
#'
#' NOTE: 00_simulate_data.R's own auto-run block is guarded by
#' `sys.nframe() == 0L || identical(environment(), globalenv())`, which is
#' true when that script is run standalone (Rscript R/00_simulate_data.R)
#' but FALSE when it is source()d from inside this function with
#' local = TRUE (as here) -- so that guard never fires in this call path.
#' source() still defines simulate_redhot() either way; we call it and save
#' explicitly rather than relying on the guarded block.
load_sim_redhot <- function(repo_root = resolve_repo_root()) {
  rds_path <- file.path(repo_root, "data", "sim_redhot.rds")
  if (!file.exists(rds_path)) {
    source(file.path(repo_root, "R", "00_simulate_data.R"), local = TRUE)
    sim_redhot <- simulate_redhot()
    dir.create(dirname(rds_path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(sim_redhot, rds_path)
    cat(sprintf(
      "Simulated %d individuals -> %s\n", nrow(sim_redhot), rds_path
    ))
  }
  readRDS(rds_path)
}

# ---------------------------------------------------------------------------
# The following four functions factor out the age-stratified
# posterior-predictive-check plotting logic that was duplicated almost
# verbatim across example28_univariate_ama1.R and example29_univariate_msp1.R.
# Behaviour is unchanged from the originals; only the code location moved.
# ---------------------------------------------------------------------------

.central_limits <- function(x, prob = 1) {
  x <- x[is.finite(x)]
  if (!length(x)) stop("No finite values available for plotting.")
  if (prob >= 1) return(range(x))
  alpha <- (1 - prob) / 2
  as.numeric(stats::quantile(x, probs = c(alpha, 1 - alpha), na.rm = TRUE))
}

.make_breaks <- function(lim, n_bins) {
  if (lim[[1L]] == lim[[2L]]) lim <- lim + c(-0.5, 0.5)
  pad <- max(diff(lim) * 1e-6, 1e-6)
  seq(lim[[1L]] - pad, lim[[2L]] + pad, length.out = n_bins + 1L)
}

.hist_density_fixed <- function(x, breaks) {
  n_bins <- length(breaks) - 1L
  out <- numeric(n_bins)
  x <- x[is.finite(x) & x >= breaks[[1L]] & x <= breaks[[length(breaks)]]]
  if (!length(x)) return(out)
  ix <- findInterval(x, breaks, rightmost.closed = TRUE)
  keep <- ix >= 1L & ix <= n_bins
  if (!any(keep)) return(out)
  out <- tabulate(ix[keep], nbins = n_bins)
  out / (sum(out) * diff(breaks))
}

.plot_age_row_uni <- function(ac, n_class, obs_y, sim_ns, sim_sp,
                               xlim, breaks, mid) {
  obs_d <- .hist_density_fixed(obs_y, breaks)
  ns_d  <- .hist_density_fixed(sim_ns, breaks)
  sp_d  <- .hist_density_fixed(sim_sp, breaks)
  ymax  <- max(obs_d, ns_d, sp_d, na.rm = TRUE) * 1.05

  graphics::plot(NA_real_, type = "n", xlim = xlim, ylim = c(0, ymax),
                 xlab = "", ylab = "Density",
                 main = sprintf("Age %s non-spatial\nn=%d", ac, n_class))
  graphics::rect(head(breaks, -1L), 0, tail(breaks, -1L), obs_d,
                 col = "grey88", border = "grey70")
  graphics::lines(mid, ns_d, col = "#2166ac", lwd = 2)
  graphics::legend("topright", legend = c("Observed", "Simulated (non-spatial)"),
                   fill = c("grey88", NA), border = c("grey70", NA),
                   lty = c(NA, 1L), col = c(NA, "#2166ac"), lwd = c(NA, 2),
                   bty = "n", cex = 0.68)

  graphics::plot(NA_real_, type = "n", xlim = xlim, ylim = c(0, ymax),
                 xlab = "", ylab = "Density",
                 main = sprintf("Age %s spatial\nn=%d", ac, n_class))
  graphics::rect(head(breaks, -1L), 0, tail(breaks, -1L), obs_d,
                 col = "grey88", border = "grey70")
  graphics::lines(mid, sp_d, col = "#b2182b", lwd = 2)
  graphics::legend("topright", legend = c("Observed", "Simulated (spatial)"),
                   fill = c("grey88", NA), border = c("grey70", NA),
                   lty = c(NA, 1L), col = c(NA, "#b2182b"), lwd = c(NA, 2),
                   bty = "n", cex = 0.68)
  invisible(NULL)
}

#' Age-stratified posterior-predictive check grid for a univariate fit.
#'
#' Simulates from the non-spatial and spatial fits within each age band and
#' compares the resulting densities against the observed histogram.
#' Reproduces the plotting logic of example28.R / example29.R exactly.
plot_univariate_ppc_grid <- function(obs_y, age, X, age_breaks,
                                      par_nonspatial, ell_id, s_hat, beta_hat,
                                      n_target = 100000L, seed0 = 1L,
                                      axis_prob = 0.99, n_bins = 30L) {
  sim_uni <- utils::getFromNamespace(".simulate_univariate_lbm_latent_t", "mlatv")
  age_class  <- cut(age, breaks = age_breaks, include.lowest = TRUE)
  age_labels <- levels(age_class)
  age_labels <- age_labels[vapply(age_labels,
    function(a) any(age_class == a, na.rm = TRUE), logical(1))]
  xlim   <- .central_limits(obs_y, axis_prob)
  breaks <- .make_breaks(xlim, n_bins)
  mid    <- (head(breaks, -1L) + tail(breaks, -1L)) / 2

  old_par <- graphics::par(mfrow = c(length(age_labels), 2L),
                           mar = c(2.6, 2.8, 2.3, 0.8), mgp = c(1.6, 0.5, 0),
                           cex.axis = 0.7, cex.lab = 0.8, cex.main = 0.82)
  on.exit(graphics::par(old_par), add = TRUE)

  for (ac in age_labels) {
    idx     <- which(age_class == ac)
    n_c     <- length(idx)
    n_rep   <- max(1L, as.integer(ceiling(n_target / n_c)))
    rep_idx <- rep(idx, times = n_rep)
    X_big   <- X[rep_idx, , drop = FALSE]

    set.seed(seed0 + match(ac, age_labels))
    sim_ns <- sim_uni(par_nonspatial, X_big)$y

    set.seed(seed0 + 1000L + match(ac, age_labels))
    eta_extra_big <- beta_hat * s_hat[ell_id[rep_idx]]
    sim_sp <- sim_uni(par_nonspatial, X_big, eta_extra = eta_extra_big)$y

    .plot_age_row_uni(ac, n_c, obs_y[idx], sim_ns, sim_sp, xlim, breaks, mid)
  }
  invisible(NULL)
}
