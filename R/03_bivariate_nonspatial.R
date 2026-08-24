# Purpose: Fit the bivariate latent-variable mixture model in a non-spatial
#   setting (Section 3 of the paper: age-dependent component locations,
#   sequential-logistic mixing probabilities, within-component correlation
#   rho_T). This corresponds to what the other scripts call "example19" -- the
#   bivariate age-delta fit that examples 22/23/28/29/31 all load as a
#   starting point. example19's own script has still not been supplied, but
#   sim_redhot_v2.R confirms the real function name and call signature
#   (fit_bivariate_lbm_latent_t_age_delta_mle_v3, with explicit breaks/mids
#   rather than a bare quadrature size M), so this version uses that instead
#   of the earlier guess (fit_bivariate_lbm_latent_t_mle).
#
# breaks/mids: built via mlatv:::.make_t_grid(M = ...), the same construction
#   example31.R uses to attach a grid to a loaded fit object.
#
# The age-spline covariate construction (Eq. age-spline in the paper) is
# reproduced from the paper text, not from example19's source, since that
# script is still not available.
#
# Output: output/rds/03_bivariate_nonspatial.rds

library(mlatv)
source(file.path(if (exists("resolve_repo_root")) resolve_repo_root() else getwd(),
                  "R", "utils_helpers.R"))
repo_root <- resolve_repo_root()
rds_dir  <- file.path(repo_root, "output", "rds")
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Configuration ----------------------------------------------------------
quadrature_M <- 12L
age_breaks   <- c(1, 5, 10, 15, 20, 40, 100)   # exploratory age bands, Sec. 5
age_knot     <- 10                              # spline knot, Eq. (age-spline)
nonspatial_control <- list(eval.max = 400L, iter.max = 400L)

# ---- Prepare data -------------------------------------------------------------
dat <- load_sim_redhot(repo_root)
dat <- dat[
  is.finite(dat$age) & dat$age > 1 &
    is.finite(dat$ama_norm) & dat$ama_norm > 0 &
    is.finite(dat$msp_norm) & dat$msp_norm > 0 &
    is.finite(dat$long_m) & is.finite(dat$lat_m),
  , drop = FALSE
]
y_1 <- log(dat$ama_norm)
y_2 <- log(dat$msp_norm)
age <- dat$age

# Linear spline in log(age) with a knot at age_knot, Eq. (age-spline) in the
# paper: d_i = (1, log(a_i), {log(a_i) - log(knot)}_+)
X <- cbind(
  `(Intercept)`   = 1,
  log_age         = log(age),
  log_age_above10 = pmax(log(age) - log(age_knot), 0)
)
delta_age <- log(age)   # covariate entering delta(a_i), Eq. (mixing_probs)

t_grid <- mlatv:::.make_t_grid(M = quadrature_M)

# ---- Fit ----------------------------------------------------------------------
fit_bivariate <- fit_bivariate_lbm_latent_t_age_delta_mle_v3(
  y_1 = y_1, y_2 = y_2, X = X, delta_age = delta_age,
  breaks = t_grid$breaks, mids = t_grid$mids,
  start = NULL, optimizer = "nlminb", control = nonspatial_control
)
if (!is.finite(fit_bivariate$logLik))
  stop("Bivariate non-spatial fit returned a non-finite log-likelihood.")
if (fit_bivariate$convergence != 0L)
  warning("Bivariate non-spatial fit did not converge: ", fit_bivariate$message)

cat("Bivariate non-spatial fit\n")
cat("logLik:", fit_bivariate$logLik, "  convergence:", fit_bivariate$convergence, "\n")
print(round(fit_bivariate$par, 4))

# ---- Save ----------------------------------------------------------------------
saveRDS(
  list(
    fit = fit_bivariate,
    data = list(y_1 = y_1, y_2 = y_2, X = X, delta_age = delta_age, age = age,
                loc_m = as.matrix(dat[, c("long_m", "lat_m")]))
  ),
  file.path(rds_dir, "03_bivariate_nonspatial.rds")
)
cat("\nSaved: output/rds/03_bivariate_nonspatial.rds\n")
