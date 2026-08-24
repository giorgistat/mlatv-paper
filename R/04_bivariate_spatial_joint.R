# Purpose: Fit the bivariate geostatistical latent variable model (Section 4),
#   both with the non-spatial parameters held fixed at the 03 fit ("Fixed
#   theta0", analogous to the paper's parameter-fixed cross-validation
#   design in Section 5.4) and fully jointly. Reports the parameter table in
#   the style of Table `tab:bootstrap_redhot`.
# Adapted from: example23_joint_spatial_z_realdata.R
# Changes from the original:
#   - data("redhot") replaced by the simulated dataset.
#   - Starting values are taken from 03_bivariate_nonspatial.rds (this
#     repo's own upstream fit) instead of the hardcoded numeric estimates
#     in example23.R, which came from the real-data fit and are meaningless
#     here.
#   - The bootstrap uncertainty step is unchanged in spirit but capped much
#     lower by default (N_BOOT) since this is meant to run on a laptop; raise
#     it for a real run.
#
# Output: output/rds/04_bivariate_spatial_joint.rds,
#         output/plot/04_spatial_z_fields.pdf

library(mlatv)
suppressPackageStartupMessages({ library(Matrix); library(fmesher) })
source(file.path(if (exists("resolve_repo_root")) resolve_repo_root() else getwd(),
                  "R", "utils_helpers.R"))
repo_root <- resolve_repo_root()
plot_dir <- file.path(repo_root, "output", "plot")
rds_dir  <- file.path(repo_root, "output", "rds")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir,  showWarnings = FALSE, recursive = TRUE)

# ---- Run settings -------------------------------------------------------------
N_SUB          <- 400L    # unique locations to keep for a fast demo run;
                           # set to NULL to use every filtered location
MAX_PER_LOC    <- 8L
M_GRID         <- 12L
MESH_CUTOFF    <- 100      # metres
MESH_MAX_EDGE  <- c(250, 700)
CTRL_FIXED     <- list(eval.max = 200, iter.max = 200)
CTRL_JOINT     <- list(eval.max = 300, iter.max = 300)
N_BOOT         <- as.integer(Sys.getenv("N_BOOT_04", "0"))   # 0 = skip
# The full-observed-information / bootstrap step is expensive; opt in with
# N_BOOT_04=200 Rscript R/04_bivariate_spatial_joint.R (or similar).

set.seed(20260824L)

# ---- 1. Load 03's fit for starting values --------------------------------------
prev_path <- file.path(rds_dir, "03_bivariate_nonspatial.rds")
if (!file.exists(prev_path)) stop("Run 03_bivariate_nonspatial.R first.")
prev <- readRDS(prev_path)
fit03 <- prev$fit

# ---- 2. Prepare data, collapse to unique locations -----------------------------
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

prep <- example_prepare_locations(dat, n_sub = N_SUB, max_per_loc = MAX_PER_LOC)
dat        <- prep$dat
loc        <- prep$loc
loc_unique <- prep$loc_unique
age        <- prep$age
y_1        <- prep$y_1
y_2        <- prep$y_2

age_knot <- 10
X <- cbind(
  `(Intercept)`   = 1,
  log_age         = log(age),
  log_age_above10 = pmax(log(age) - log(age_knot), 0)
)
delta_age <- log(age)

# ---- 3. Build mesh --------------------------------------------------------------
mesh <- fmesher::fm_mesh_2d_inla(loc = loc_unique, max.edge = MESH_MAX_EDGE,
                                  cutoff = MESH_CUTOFF)
cat("Mesh nodes:", mesh$n, "\n")
example_save_mesh_plot(mesh, dat, loc_unique, site_center,
                        prep$n_unique, prep$n_obs,
                        file.path(plot_dir, "04_mesh.pdf"))

# ---- 4. Fixed-theta0 fit (spatial parameters only, theta0 held at 03) ---------
theta0_start <- as.list(fit03$par)
spatial_start <- list(range_1 = 300, range_2 = 300, rho = 0.3,
                       beta_A = 1, beta_M = 1)
start_full <- c(list(theta0 = theta0_start), spatial_start)

cat("\n--- Fixed-theta0 fit (estimate_theta0 = FALSE) ---\n")
t0 <- proc.time()
fit_fixed <- fit_bivariate_spatial_z_joint(
  y_1 = y_1, y_2 = y_2, X = X, delta_age = delta_age,
  mesh = mesh, loc = loc, M = M_GRID,
  start = start_full, estimate_theta0 = FALSE,
  control = CTRL_FIXED, verbose = 1L
)
t_fixed <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("Fixed-theta0: logLik=%.2f  convergence=%d  time=%.1fs\n",
            fit_fixed$logLik, fit_fixed$convergence, t_fixed))

# ---- 5. Joint fit (all parameters jointly), warm-started from step 4 -----------
start_joint <- c(
  list(theta0 = theta0_start),
  as.list(fit_fixed$par[c("range_1", "range_2", "rho", "beta_A", "beta_M")])
)
cat("\n--- Joint fit (estimate_theta0 = TRUE) ---\n")
t0 <- proc.time()
fit_joint <- fit_bivariate_spatial_z_joint(
  y_1 = y_1, y_2 = y_2, X = X, delta_age = delta_age,
  mesh = mesh, loc = loc, M = M_GRID,
  start = start_joint, estimate_theta0 = TRUE,
  control = CTRL_JOINT, verbose = 1L
)
t_joint <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("Joint fit:    logLik=%.2f  convergence=%d  time=%.1fs\n",
            fit_joint$logLik, fit_joint$convergence, t_joint))

# ---- 6. Parameter table (style of Table tab:bootstrap_redhot) -----------------
par_table <- data.frame(
  parameter = names(fit_joint$par),
  estimate  = as.numeric(fit_joint$par),
  row.names = NULL
)
cat("\nFitted spatial parameters (joint fit):\n")
print(round(par_table$estimate[match(
  c("range_1", "range_2", "rho", "beta_A", "beta_M"), par_table$parameter
)], 3))

# ---- 7. Bootstrap uncertainty (opt-in, expensive) ------------------------------
uncertainty_joint <- NULL
if (N_BOOT > 0L) {
  cat(sprintf("\nComputing observed information (N_BOOT=%d)...\n", N_BOOT))
  uncertainty_joint <- tryCatch(
    parameter_uncertainty(fit_joint, level = 0.95, scope = "joint",
                           control = list(progress = TRUE)),
    error = function(e) {
      warning("Uncertainty calculation failed: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(uncertainty_joint)) print(uncertainty_joint$summary, row.names = FALSE)
} else {
  cat("\nSkipping bootstrap uncertainty (set N_BOOT_04 > 0 to enable).\n")
}

# ---- 8. Posterior mean field plots ---------------------------------------------
s_fixed <- fit_fixed$s_hat
s_joint <- fit_joint$s_hat
grDevices::pdf(file.path(plot_dir, "04_spatial_z_fields.pdf"), width = 10, height = 4)
graphics::par(mfrow = c(1, 4), mar = c(3, 3, 2, 1))
for (nm in c("S_A", "S_M")) {
  for (flab in c("fixed", "joint")) {
    sv <- if (flab == "fixed") s_fixed[, nm] else s_joint[, nm]
    clim <- max(abs(sv)); clim <- c(-clim, clim)
    graphics::plot(loc_unique,
      col = grDevices::colorRampPalette(c("blue", "white", "red"))(100)[
        pmin(100, pmax(1, round((sv - clim[1]) / diff(clim) * 99) + 1))],
      pch = 16, cex = 0.4,
      main = sprintf("%s (%s)", nm, flab), xlab = "x (m)", ylab = "y (m)")
  }
}
graphics::par(mfrow = c(1, 1))
grDevices::dev.off()
cat("Saved:", file.path(plot_dir, "04_spatial_z_fields.pdf"), "\n")

# ---- 9. Save ---------------------------------------------------------------------
saveRDS(
  list(
    fit_fixed = fit_fixed, fit_joint = fit_joint,
    uncertainty_joint = uncertainty_joint,
    par_table = par_table,
    mesh = mesh, loc = loc, loc_unique = loc_unique,
    X = X, y_1 = y_1, y_2 = y_2, age = age, delta_age = delta_age,
    site_center = site_center,
    t_fixed = t_fixed, t_joint = t_joint
  ),
  file.path(rds_dir, "04_bivariate_spatial_joint.rds")
)
cat("Saved: output/rds/04_bivariate_spatial_joint.rds\n")
