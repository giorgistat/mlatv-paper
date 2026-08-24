###############################################################################
## 07_simulation_study.R
##
## Adapted from the REAL sim_redhot_v2.R (supplied after this repo's first
## draft, which only had a documented stub here). The loop logic, formulas,
## grid, and output schema below are ported essentially unchanged from the
## original; the only substantive change is where the "real domain" and
## "truth" fit come from:
##
##   ORIGINAL                              THIS REPO
##   readRDS("article_data.rds")$spatial   built from 04's saved mesh/X/age/
##                                          delta_age/loc/loc_unique (which
##                                          in turn come from the simulated
##                                          dataset, not the real survey)
##   readRDS("spatial_joint_fit.rds")      04_bivariate_spatial_joint.rds's
##                                          fit_joint object
##
## Everything downstream of that substitution -- the N_SAMPLE x rho_S x
## rho_T x sigma_scale loop, the Laplace re-draw of the truth field, the
## three fitting arms, the S/T/Y prediction step, and the saved schema -- is
## the original design, not a stub.
##
## Design: N_SAMPLE x rho_S x rho_T x sigma_scale, with N_SAMPLE in
## {300, 600, 1200} and rho_S, rho_T each in {0.1, 0.5, 0.9}, using the
## fitted ranges (range_1, range_2) from the reference fit -- rho_S is
## therefore the coupling coefficient, not necessarily the pointwise
## correlation, since the two ranges need not be equal.
##
## NO SEEDING (preserved from the original): every replicate draws a fresh
## random subsample of locations/residents, a fresh spatial field, and fresh
## simulated responses -- nothing is fixed via set.seed() anywhere in the
## main loop, so repeated/parallel runs produce independent replicates
## without needing to coordinate seeds across tasks.
##
## DEMO_MODE (added for this repo only, not in the original): the full grid
## below is 3 x 3 x 3 x 4 = 108 scenarios, each fitting three spatial models
## -- far too slow for a laptop smoke test. Set DEMO_MODE = FALSE to run the
## real published grid (e.g. as one task of a SLURM array job); the default
## here (TRUE) shrinks every grid to length 1 so the script exercises the
## full pipeline once, quickly.
##
## Output per run: true and predicted S, T, Y at every sampled location for
## all three arms (with N_SAMPLE recorded on every row), plus a one-row
## parameter-estimate summary per replicate (also carrying N_SAMPLE). Saved
## as one .rds; the random integer used in the file name is generated at the
## very end of the script, after all computation, purely for naming.
###############################################################################
library(mlatv)
library(Matrix)
source(file.path(if (exists("resolve_repo_root")) resolve_repo_root() else getwd(),
                  "R", "utils_helpers.R"))
repo_root <- resolve_repo_root()
rds_dir  <- file.path(repo_root, "output", "rds")
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

DEMO_MODE <- tolower(Sys.getenv("SIM07_FULL_GRID", "false")) != "true"
# DEMO_MODE is TRUE (shrunk grid) unless SIM07_FULL_GRID=true is set.

## ---- run settings ---------------------------------------------------------
N_SAMPLE_GRID <- if (DEMO_MODE) c(300L) else c(300L, 600L, 1200L)
## Number of individuals simulated per sampled location (household/cluster
## size). N_SAMPLE is the TOTAL number of individuals; the number of distinct
## locations sampled is N_SAMPLE %/% N_PER_LOC. Requires N_SAMPLE to be an
## exact multiple of N_PER_LOC -- checked below.
N_PER_LOC     <- 5L
RHO_S_GRID    <- if (DEMO_MODE) c(0.5) else c(0.1, 0.5, 0.9)
RHO_T_GRID    <- if (DEMO_MODE) c(0.5) else c(0.1, 0.5, 0.9)
## Scales the Y | T observation-model VARIANCE: sigma^2_sim = SIGMA_SCALE *
## sigma^2_fitted for all four conditional SDs (sigma0_1, sigma1_1, sigma0_2,
## sigma1_2). 1 reproduces the fitted noise level exactly, as a check that
## this addition is a no-op there; smaller values raise the signal-to-noise
## ratio of Y as a proxy for T, without touching the latent process itself.
SIGMA_SCALE_GRID <- if (DEMO_MODE) c(1) else c(1, 0.5, 0.25, 0.1)
N_REP         <- 1L

if (any(N_SAMPLE_GRID %% N_PER_LOC != 0L))
  stop("Every value in N_SAMPLE_GRID must be an exact multiple of N_PER_LOC (",
       N_PER_LOC, "). N_SAMPLE_GRID = ", paste(N_SAMPLE_GRID, collapse = ", "),
       call. = FALSE)

MLE_CONTROL     <- list(eval.max = 2000L, iter.max = 1000L)
SPATIAL_CONTROL <- list(eval.max = 400L, iter.max = 400L)
LAPLACE_CONTROL <- list(tol = 1e-5)
VERBOSE <- 0L

## ---- build the domain and reference ("truth") fit from THIS repo's own
## upstream fit, instead of article_data.rds / spatial_joint_fit.rds --------
fit04_path <- file.path(rds_dir, "04_bivariate_spatial_joint.rds")
if (!file.exists(fit04_path))
  stop("Run 04_bivariate_spatial_joint.R first (its fit_joint is used as ",
       "the simulation-study 'truth' fit, in place of the original's ",
       "spatial_joint_fit.rds).")
fit04 <- readRDS(fit04_path)
truth_fit <- fit04$fit_joint

age_breaks <- c(1, 5, 10, 15, 20, 40, 100)
age_class_full <- cut(fit04$age, breaks = age_breaks, include.lowest = TRUE)
t_grid <- mlatv:::.make_t_grid(M = 12L)   # rebuilt explicitly rather than
                                           # trusting inherited breaks/mids,
                                           # matching example31.R's pattern

spatial_data <- list(
  mesh       = fit04$mesh,
  loc_unique = fit04$loc_unique,
  age        = fit04$age,
  age_class  = age_class_full,
  X          = fit04$X,
  delta_age  = fit04$delta_age,
  n_unique   = nrow(fit04$loc_unique)
)

spatial_parameter_names <- c("range_1", "range_2", "rho", "beta_A", "beta_M")
range_1_true <- unname(truth_fit$par[["range_1"]])
range_2_true <- unname(truth_fit$par[["range_2"]])
beta_A_true  <- unname(truth_fit$par[["beta_A"]])
beta_M_true  <- unname(truth_fit$par[["beta_M"]])
cat(sprintf("Reference fit: range_1=%.2f  range_2=%.2f  beta_A=%.3f  beta_M=%.3f\n",
            range_1_true, range_2_true, beta_A_true, beta_M_true))
cat(sprintf("DEMO_MODE=%s -> grid is %d x %d x %d x %d = %d scenarios\n",
            DEMO_MODE, length(N_SAMPLE_GRID), length(RHO_S_GRID),
            length(RHO_T_GRID), length(SIGMA_SCALE_GRID),
            length(N_SAMPLE_GRID) * length(RHO_S_GRID) *
              length(RHO_T_GRID) * length(SIGMA_SCALE_GRID)))

## =================================================================================
## Draw one Laplace-posterior sample of the FEM weights, w | y ~ N(w_hat, H_hat^{-1}).
## =================================================================================
sample_laplace_w <- function(w_hat, H) {
  H  <- as(Matrix::forceSymmetric(H), "CsparseMatrix")
  ch <- Matrix::Cholesky(H, perm = TRUE, LDL = FALSE)
  d  <- length(w_hat)
  z  <- stats::rnorm(d)
  pert <- Matrix::solve(ch, Matrix::solve(ch, z, system = "Lt"), system = "Pt")
  w_hat + as.numeric(pert)
}

## =================================================================================
## Predict S, T, Y at the sampled locations from the fitted JOINT bivariate model:
## one Laplace-posterior draw of (S1, S2) -> forward-simulate via the package's
## own .simulate_bivariate_lbm_latent_t().
## =================================================================================
predict_joint_STY <- function(fit_spatial, X, delta_age) {
  bsz   <- fit_spatial$bsz_data
  n_loc <- bsz$n_loc_unique
  w_draw <- sample_laplace_w(fit_spatial$w_hat, fit_spatial$H_hat)
  s_loc  <- as.numeric(bsz$A_block %*% w_draw)
  S1_loc <- s_loc[seq_len(n_loc)]; S2_loc <- s_loc[n_loc + seq_len(n_loc)]
  S1_i <- S1_loc[bsz$ell_id]; S2_i <- S2_loc[bsz$ell_id]

  beta_A <- unname(fit_spatial$par[["beta_A"]])
  beta_M <- unname(fit_spatial$par[["beta_M"]])

  sim <- mlatv:::.simulate_bivariate_lbm_latent_t(
    par = fit_spatial$theta0_par, X = X, delta_age = delta_age,
    eta_extra_1 = beta_A * S1_i, eta_extra_2 = beta_M * S2_i,
    breaks = fit_spatial$breaks, mids = fit_spatial$mids,
    latent_method = "rejection"
  )
  data.frame(S1 = S1_i, S2 = S2_i, T1 = sim$T_1, T2 = sim$T_2,
             Y1 = sim$y_1, Y2 = sim$y_2)
}

## =================================================================================
## Predict S, T, Y at the sampled locations from ONE fitted univariate
## (single-antigen) spatial model: Laplace-draw S -> forward-simulate via
## .simulate_univariate_lbm_latent_t().
## =================================================================================
predict_uni_STY <- function(fit_uni, X) {
  usz   <- fit_uni$usz_data
  w_draw <- sample_laplace_w(fit_uni$w_hat, fit_uni$H_hat)
  s_loc  <- as.numeric(usz$A %*% w_draw)
  S_i    <- s_loc[usz$ell_id]

  beta_k <- unname(fit_uni$par[["beta"]])
  sim <- mlatv:::.simulate_univariate_lbm_latent_t(
    par = fit_uni$nonspatial_fit$par, X = X, eta_extra = beta_k * S_i
  )
  data.frame(S = S_i, T = sim$T, Y = sim$y)
}

## =================================================================================
## Main N_SAMPLE x scenario x replicate loop (unchanged from the original)
## =================================================================================
records     <- list()
predictions <- list()
run_i <- 0L

for (N_SAMPLE in N_SAMPLE_GRID) {
  for (rho_S in RHO_S_GRID) {
    for (rho_T in RHO_T_GRID) {
      for (sigma_scale in SIGMA_SCALE_GRID) {
        for (rep_i in seq_len(N_REP)) {
          run_i <- run_i + 1L
          cat(sprintf("\n=== n_sample=%d rho_S=%.1f rho_T=%.1f sigma_scale=%.2f rep=%d ===\n",
                      N_SAMPLE, rho_S, rho_T, sigma_scale, rep_i))

          ## ---- (0) fresh random subsample of locations, N_PER_LOC individuals
          ## per location -----------------------------------------------------
          n_loc <- N_SAMPLE %/% N_PER_LOC
          location_id_unique <- sample.int(spatial_data$n_unique, size = n_loc)

          row_id <- unlist(lapply(location_id_unique, function(id) {
            candidates <- which(truth_fit$bsz_data$ell_id == id)
            sample(candidates, size = N_PER_LOC,
                   replace = length(candidates) < N_PER_LOC)
          }))
          ## individual-level location index, repeated N_PER_LOC times per
          ## sampled location, in the same order as row_id
          location_id <- rep(location_id_unique, each = N_PER_LOC)

          loc       <- spatial_data$loc_unique[location_id, , drop = FALSE]
          age       <- spatial_data$age[row_id]
          age_class <- spatial_data$age_class[row_id]
          X         <- spatial_data$X[row_id, , drop = FALSE]
          delta_age <- spatial_data$delta_age[row_id]

          ## ---- (1) truth for this scenario: reference theta0, scenario-specific
          ## rho_T/rho_S -------------------------------------------------------
          theta0_truth <- truth_fit$theta0_par
          theta0_truth[["rho_T"]] <- rho_T
          ## Scale the four Y | T conditional SDs by sqrt(sigma_scale), so that
          ## the conditional VARIANCE scales by sigma_scale exactly; the
          ## component means mu0/mu1 and the latent-T structure are untouched.
          sigma_fields <- c("sigma0_1", "sigma1_1", "sigma0_2", "sigma1_2")
          missing_sigma <- setdiff(sigma_fields, names(theta0_truth))
          if (length(missing_sigma)) {
            stop("theta0_par is missing expected sigma fields: ",
                 paste(missing_sigma, collapse = ", "),
                 ". Check names(truth_fit$theta0_par) and update sigma_fields.",
                 call. = FALSE)
          }
          for (nm in sigma_fields) {
            theta0_truth[[nm]] <- theta0_truth[[nm]] * sqrt(sigma_scale)
          }
          theta0_start <- mlatv:::.bivariate_lbm_latent_t_age_delta_vector_to_list(
            truth_fit$theta0_par_working, X
          )

          spatial_truth <- truth_fit$par[spatial_parameter_names]
          spatial_truth[["rho"]] <- rho_S

          ## ---- (2) simulate the spatial fields over the FULL domain, then subset
          field <- mlatv::simulate_bivariate_field(
            mesh = spatial_data$mesh, loc = spatial_data$loc_unique,
            range1 = spatial_truth[["range_1"]], range2 = spatial_truth[["range_2"]],
            sigma1 = 1, sigma2 = 1, rho_S = spatial_truth[["rho"]]
          )

          ## ---- (3) simulate one paired response at each sampled location -------
          simulated <- mlatv:::.simulate_bivariate_lbm_latent_t(
            par = theta0_truth, X = X, delta_age = delta_age,
            eta_extra_1 = spatial_truth[["beta_A"]] * field$S1[location_id],
            eta_extra_2 = spatial_truth[["beta_M"]] * field$S2[location_id],
            breaks = t_grid$breaks, mids = t_grid$mids,
            latent_method = "rejection"
          )
          y_1 <- simulated$y_1; y_2 <- simulated$y_2

          ## ---- (4) ARM 1: bivariate joint model, two-stage full MLE -------------
          fit_biv <- tryCatch({
            fit_mle <- fit_bivariate_lbm_latent_t_age_delta_mle_v3(
              y_1 = y_1, y_2 = y_2, X = X, delta_age = delta_age,
              breaks = t_grid$breaks, mids = t_grid$mids,
              start = theta0_start, optimizer = "nlminb", control = MLE_CONTROL
            )
            fit_spatial <- fit_bivariate_spatial_z_joint(
              y_1 = y_1, y_2 = y_2, X = X, delta_age = delta_age,
              mesh = spatial_data$mesh, loc = loc,
              breaks = fit_mle$breaks, mids = fit_mle$mids,
              start = as.list(spatial_truth), estimate_theta0 = FALSE,
              nonspatial_fit = fit_mle, nonspatial_convergence = "warn",
              field_structure = "correlated",
              control = SPATIAL_CONTROL, laplace_control = LAPLACE_CONTROL,
              information = "none", verbose = VERBOSE
            )
            list(mle = fit_mle, spatial = fit_spatial)
          }, error = function(e) { message("bivariate fit failed: ", conditionMessage(e)); NULL })

          ## ---- (5) ARM 2 & 3: univariate spatial model per antigen, two-stage full MLE
          fit_uni1 <- tryCatch(
            fit_univariate_lbm_spatial(
              y = y_1, mesh = spatial_data$mesh, loc = loc, X = X,
              nonspatial_method = "mle", nonspatial_convergence = "warn",
              breaks = t_grid$breaks, mids = t_grid$mids,
              spatial_control = SPATIAL_CONTROL, laplace_control = LAPLACE_CONTROL,
              information = "none", verbose = VERBOSE
            ),
            error = function(e) { message("univariate antigen-1 fit failed: ",
                                          conditionMessage(e)); NULL })
          fit_uni2 <- tryCatch(
            fit_univariate_lbm_spatial(
              y = y_2, mesh = spatial_data$mesh, loc = loc, X = X,
              nonspatial_method = "mle", nonspatial_convergence = "warn",
              breaks = t_grid$breaks, mids = t_grid$mids,
              spatial_control = SPATIAL_CONTROL, laplace_control = LAPLACE_CONTROL,
              information = "none", verbose = VERBOSE
            ),
            error = function(e) { message("univariate antigen-2 fit failed: ",
                                          conditionMessage(e)); NULL })

          ## ---- (6) predict S, T, Y at the sampled locations for all three arms ----
          pj <- if (!is.null(fit_biv))
            tryCatch(predict_joint_STY(fit_biv$spatial, X, delta_age),
                     error = function(e) { message("joint prediction failed: ",
                                                   conditionMessage(e)); NULL })
          p1 <- if (!is.null(fit_uni1))
            tryCatch(predict_uni_STY(fit_uni1, X),
                     error = function(e) { message("uni1 prediction failed: ",
                                                   conditionMessage(e)); NULL })
          p2 <- if (!is.null(fit_uni2))
            tryCatch(predict_uni_STY(fit_uni2, X),
                     error = function(e) { message("uni2 prediction failed: ",
                                                   conditionMessage(e)); NULL })

          na_col <- rep(NA_real_, N_SAMPLE)
          pred_table <- data.frame(
            run = run_i, n_sample = N_SAMPLE,
            rho_S_true = rho_S, rho_T_true = rho_T, sigma_scale = sigma_scale,
            rep = rep_i,
            location_id = location_id, source_row = row_id,
            loc_x = loc[, 1L], loc_y = loc[, 2L], age = age,
            ## true
            S1_true = field$S1[location_id], S2_true = field$S2[location_id],
            T1_true = simulated$T_1, T2_true = simulated$T_2,
            Y1_true = y_1, Y2_true = y_2,
            ## bivariate joint predicted
            S1_joint = if (!is.null(pj)) pj$S1 else na_col,
            S2_joint = if (!is.null(pj)) pj$S2 else na_col,
            T1_joint = if (!is.null(pj)) pj$T1 else na_col,
            T2_joint = if (!is.null(pj)) pj$T2 else na_col,
            Y1_joint = if (!is.null(pj)) pj$Y1 else na_col,
            Y2_joint = if (!is.null(pj)) pj$Y2 else na_col,
            ## separate univariate predicted
            S1_sep = if (!is.null(p1)) p1$S else na_col,
            T1_sep = if (!is.null(p1)) p1$T else na_col,
            Y1_sep = if (!is.null(p1)) p1$Y else na_col,
            S2_sep = if (!is.null(p2)) p2$S else na_col,
            T2_sep = if (!is.null(p2)) p2$T else na_col,
            Y2_sep = if (!is.null(p2)) p2$Y else na_col
          )
          predictions[[run_i]] <- pred_table

          ## ---- (7) collect parameter estimates ----------------------------------
          grab <- function(fit, nm) if (!is.null(fit) && nm %in% names(fit$par)) fit$par[[nm]] else NA_real_
          rec <- data.frame(
            run = run_i, n_sample = N_SAMPLE,
            rho_S_true = rho_S, rho_T_true = rho_T, sigma_scale = sigma_scale,
            rep = rep_i,
            biv_range_1 = if (!is.null(fit_biv)) grab(fit_biv$spatial, "range_1") else NA_real_,
            biv_range_2 = if (!is.null(fit_biv)) grab(fit_biv$spatial, "range_2") else NA_real_,
            biv_rho     = if (!is.null(fit_biv)) grab(fit_biv$spatial, "rho")     else NA_real_,
            biv_beta_A  = if (!is.null(fit_biv)) grab(fit_biv$spatial, "beta_A") else NA_real_,
            biv_beta_M  = if (!is.null(fit_biv)) grab(fit_biv$spatial, "beta_M") else NA_real_,
            biv_rho_T   = if (!is.null(fit_biv)) fit_biv$spatial$theta0_par[["rho_T"]] else NA_real_,
            biv_logLik  = if (!is.null(fit_biv)) fit_biv$spatial$logLik else NA_real_,
            biv_conv    = if (!is.null(fit_biv)) fit_biv$spatial$convergence else NA_integer_,
            biv_mle_conv = if (!is.null(fit_biv)) fit_biv$mle$convergence else NA_integer_,
            uni1_range = if (!is.null(fit_uni1)) grab(fit_uni1, "range") else NA_real_,
            uni1_beta  = if (!is.null(fit_uni1)) grab(fit_uni1, "beta")  else NA_real_,
            uni1_conv  = if (!is.null(fit_uni1)) fit_uni1$convergence else NA_integer_,
            uni2_range = if (!is.null(fit_uni2)) grab(fit_uni2, "range") else NA_real_,
            uni2_beta  = if (!is.null(fit_uni2)) grab(fit_uni2, "beta")  else NA_real_,
            uni2_conv  = if (!is.null(fit_uni2)) fit_uni2$convergence else NA_integer_
          )
          records[[run_i]] <- rec

          cat(sprintf("  biv conv=%s (mle=%s)  uni1 conv=%s  uni2 conv=%s\n",
                      rec$biv_conv, rec$biv_mle_conv, rec$uni1_conv, rec$uni2_conv))
        }
      }
    }
  }
}

estimates_df   <- do.call(rbind, records)
predictions_df <- do.call(rbind, predictions)

## ---- save: the random integer used in the file name is generated here, at
## the very end, purely for a collision-free file name -- it is not a seed
## and has no bearing on anything computed above.
run_tag  <- sample.int(.Machine$integer.max, 1L)
out_file <- file.path(rds_dir, sprintf("07_simulation_study_%d.rds", run_tag))

saveRDS(
  list(
    run_tag     = run_tag,
    demo_mode   = DEMO_MODE,
    settings    = list(N_SAMPLE_GRID = N_SAMPLE_GRID,
                       RHO_S_GRID = RHO_S_GRID, RHO_T_GRID = RHO_T_GRID,
                       SIGMA_SCALE_GRID = SIGMA_SCALE_GRID,
                       sigma_fields = c("sigma0_1", "sigma1_1",
                                        "sigma0_2", "sigma1_2"),
                       sigma_fitted = truth_fit$theta0_par[
                         c("sigma0_1", "sigma1_1", "sigma0_2", "sigma1_2")],
                       N_REP = N_REP,
                       range_1_true = range_1_true, range_2_true = range_2_true,
                       beta_A_true = beta_A_true, beta_M_true = beta_M_true,
                       reference_fit_source = fit04_path),
    estimates   = estimates_df,     # one row per run (has n_sample column)
    predictions = predictions_df    # one row per sampled location per run
    # (has n_sample column) -- ready to rbind across parallel-run files
  ),
  out_file
)
cat(sprintf("\nSaved: %s\n", out_file))
cat("Done.\n")
cat(
  "\nNOTE: this ran with DEMO_MODE=", DEMO_MODE,
  ". Set the environment variable SIM07_FULL_GRID=true before running this ",
  "script to use the full published grid (108 scenarios); expect this to ",
  "take a long time and to be best run as a SLURM/HPC array job, one task ",
  "per (scenario, replicate), as in the original.\n", sep = ""
)
