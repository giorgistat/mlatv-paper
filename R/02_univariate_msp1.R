# Purpose: Fit the univariate latent-T model and its spatial-Z extension to
#   MSP1 (Section 5.1 of the paper). Adapted from example28_univariate_ama1.R.
# Changes from the original:
#   - data("redhot") replaced by the simulated dataset (load_sim_redhot()).
#   - The hardcoded `nonspatial_start` (natural-scale values from a prior fit
#     on the real data) is dropped; the fit starts from package defaults
#     (start = NULL). This will be slower to converge and is NOT expected to
#     reproduce the paper's MSP1 estimates -- only to exercise the pipeline.
#     When re-running on the real data, restore your own start values.
# Output: output/rds/02_univariate_msp1.rds, output/plot/02_msp1_age_histograms.pdf

library(mlatv)

source_root <- if (exists("resolve_repo_root")) resolve_repo_root() else getwd()
source(file.path(source_root, "R", "utils_helpers.R"))
repo_root <- resolve_repo_root()
plot_dir <- file.path(repo_root, "output", "plot")
rds_dir  <- file.path(repo_root, "output", "rds")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir,  showWarnings = FALSE, recursive = TRUE)

# ---- Configuration ------------------------------------------------------
quadrature_M  <- 8L
age_breaks    <- c(1, 5, 10, 15, 20, 40, 100)
mesh_range    <- 650
mesh_cutoff   <- 180
mesh_max_edge <- c(260, 650)
verbose       <- 1L
nonspatial_control <- list(eval.max = 550L, iter.max = 550L,
                            trace = if (verbose >= 1L) 1L else 0L)
spatial_control    <- list(eval.max = 40L, iter.max = 40L)

# ---- Prepare data ---------------------------------------------------------
dat <- load_sim_redhot(repo_root)
analysis_data <- dat[
  is.finite(dat$age) & dat$age > 1 &
    is.finite(dat$msp_norm) & dat$msp_norm > 0 &
    is.finite(dat$long_m) & is.finite(dat$lat_m),
  , drop = FALSE
]

y   <- log(analysis_data$msp_norm)
age <- analysis_data$age
age_class <- cut(age, breaks = age_breaks, include.lowest = TRUE)
if (anyNA(age_class)) stop("All retained ages must fall within age_breaks.")
X <- stats::model.matrix(~ age_class)
colnames(X)[1L] <- "(Intercept)"
colnames(X)[-1L] <- paste0("age_", make.names(levels(age_class)[-1L]))
loc <- as.matrix(analysis_data[, c("long_m", "lat_m")])
storage.mode(loc) <- "double"
loc <- sweep(loc, 2L, colMeans(loc), "-")

# ---- Fit non-spatial latent-T model ---------------------------------------
fit_nonspatial <- fit_univariate_lbm_latent_t_mle(
  y = y, X = X, M = quadrature_M,
  start = NULL,               # see header note: original used a hot start
  control = nonspatial_control
)
if (!is.finite(fit_nonspatial$logLik))
  stop("MSP1 non-spatial fit returned a non-finite log-likelihood.")
if (fit_nonspatial$convergence != 0L)
  warning("MSP1 non-spatial fit did not converge: ", fit_nonspatial$message)

component_evidence <- univariate_lbm_component_evidence(fit_nonspatial)

# ---- Spatial-Z extension ---------------------------------------------------
mesh <- build_rspde_mesh(loc = loc, range = mesh_range, cutoff = mesh_cutoff,
                          max.edge = mesh_max_edge)
spatial_data <- build_univariate_spatial_z_data(
  mesh = mesh, loc = loc,
  log_q = component_evidence$log_q, eta = fit_nonspatial$eta
)
fit_spatial <- fit_univariate_spatial_z(
  usz_data = spatial_data, control = spatial_control, verbose = verbose
)
if (!is.finite(fit_spatial$logLik))
  stop("MSP1 spatial fit returned a non-finite log-likelihood.")
if (fit_spatial$convergence != 0L)
  warning("MSP1 spatial fit did not converge: ", fit_spatial$message)

# ---- Report -----------------------------------------------------------------
fit_dimensions <- data.frame(
  response = "MSP1", observations = nrow(analysis_data),
  unique_sites = spatial_data$n_loc_unique, mesh_nodes = spatial_data$n_mesh,
  quadrature_cells = length(fit_nonspatial$mids)
)
fit_status <- data.frame(
  stage = c("non-spatial", "spatial"),
  convergence = c(fit_nonspatial$convergence, fit_spatial$convergence),
  logLik = c(fit_nonspatial$logLik, fit_spatial$logLik)
)
cat("MSP1 fit dimensions\n"); print(fit_dimensions, row.names = FALSE)
cat("\nMSP1 fit status\n");   print(fit_status, row.names = FALSE)
cat("\nMSP1 observation and latent parameters\n")
print(round(fit_nonspatial$par[c("mu0", "mu1", "sigma0", "sigma1", "nu_T")], 4))
cat("\nMSP1 spatial parameters\n")
print(round(fit_spatial$par[c("range", "beta")], 4))

# ---- Age-stratified posterior-predictive check ------------------------------
fig_hist <- file.path(plot_dir, "02_msp1_age_histograms.pdf")
grDevices::pdf(fig_hist, width = 9, height = 3.1 * length(levels(age_class)))
plot_univariate_ppc_grid(
  obs_y = y, age = age, X = X, age_breaks = age_breaks,
  par_nonspatial = fit_nonspatial$par,
  ell_id = spatial_data$ell_id, s_hat = fit_spatial$s_hat,
  beta_hat = fit_spatial$par[["beta"]], seed0 = 20260102L
)
grDevices::dev.off()
cat("Saved:", fig_hist, "\n")

# ---- Save -------------------------------------------------------------------
saveRDS(
  list(
    fit_nonspatial = fit_nonspatial,
    fit_spatial    = fit_spatial,
    spatial_data   = spatial_data,
    fit_dimensions = fit_dimensions,
    fit_status     = fit_status,
    X = X, y = y, age = age, loc = loc
  ),
  file.path(rds_dir, "02_univariate_msp1.rds")
)
cat("\nSaved: output/rds/02_univariate_msp1.rds\n")
