# Purpose: Simulate a bivariate antibody dataset (sim_redhot) with the same
#   column structure, scale, and rough age/spatial dependence as the real
#   RedHot serosurvey used in the paper, so the rest of the pipeline can run
#   end-to-end without the (non-redistributable) real data.
# Output: data/sim_redhot.rds and the same object returned invisibly.
#
# NOTE: This is a structural stand-in, not a faithful draw from the fitted
#   bivariate latent-variable model. It uses a simple age-trend + spatial
#   Gaussian-process + bivariate-normal-noise construction on the log scale,
#   calibrated only to be roughly the right order of magnitude and roughly
#   right shape (positive age trend, positive spatial clustering, positive
#   cross-antigen correlation). Do not use it to check numerical agreement
#   with the paper's reported estimates -- only to exercise the pipeline.

simulate_redhot <- function(
    n_individual   = 4000L,
    n_location     = 800L,
    age_max        = 80,
    domain_km      = 6,
    spatial_range_km = 0.5,
    cross_field_rho   = 0.6,
    within_host_rho   = 0.6,
    seed = 20260824L
) {
  set.seed(seed)

  # --- locations -----------------------------------------------------------
  loc_km <- cbind(
    x = stats::runif(n_location, -domain_km / 2, domain_km / 2),
    y = stats::runif(n_location, -domain_km / 2, domain_km / 2)
  )
  dist_km <- as.matrix(stats::dist(loc_km))
  corr <- exp(-dist_km / spatial_range_km)

  chol_joint <- function(rho) {
    sigma <- rbind(
      cbind(corr, rho * corr),
      cbind(rho * corr, corr)
    )
    sigma <- sigma + diag(1e-8, nrow(sigma))
    chol(sigma)
  }
  root <- chol_joint(cross_field_rho)
  z <- as.numeric(matrix(stats::rnorm(2L * n_location), nrow = 1L) %*% root)
  s1 <- z[seq_len(n_location)]
  s2 <- z[n_location + seq_len(n_location)]
  s1 <- s1 / stats::sd(s1)
  s2 <- s2 / stats::sd(s2)

  # --- individuals: assign to locations, simulate age ----------------------
  loc_id <- sample.int(n_location, n_individual, replace = TRUE)
  age <- stats::rgamma(n_individual, shape = 1.4, rate = 1.4 / 12)
  age <- pmin(pmax(age, 1), age_max)

  # --- age-dependent mean seroreactivity on logit scale ---------------------
  logit_mean_1 <- -2.0 + 0.9 * log(age)
  logit_mean_2 <- -1.6 + 0.7 * log(age)
  spatial_effect_1 <- 1.1 * s1[loc_id]
  spatial_effect_2 <- 0.9 * s2[loc_id]

  eta1 <- logit_mean_1 + spatial_effect_1
  eta2 <- logit_mean_2 + spatial_effect_2

  # within-host correlated noise on the latent scale, then map through a
  # logistic link into (0,1) "seroreactivity", then into log antibody units
  host_noise <- MASS_mvrnorm_fallback(
    n_individual,
    rho = within_host_rho
  )
  t1 <- stats::plogis(eta1 + 0.6 * host_noise[, 1L])
  t2 <- stats::plogis(eta2 + 0.6 * host_noise[, 2L])

  mu0_1 <- -6.5; mu1_1 <- 0.9; sigma0_1 <- 0.55; sigma1_1 <- 0.20
  mu0_2 <- -6.9; mu1_2 <- 0.8; sigma0_2 <- 0.60; sigma1_2 <- 0.30

  mean1 <- mu0_1 + t1 * (mu1_1 - mu0_1)
  mean2 <- mu0_2 + t2 * (mu1_2 - mu0_2)
  sd1   <- sqrt((1 - t1) * sigma0_1^2 + t1 * sigma1_1^2)
  sd2   <- sqrt((1 - t2) * sigma0_2^2 + t2 * sigma1_2^2)

  log_y1 <- mean1 + stats::rnorm(n_individual, 0, sd1)
  log_y2 <- mean2 + stats::rnorm(n_individual, 0, sd2)

  ama_norm <- exp(log_y1)
  msp_norm <- exp(log_y2)

  loc_m <- loc_km[loc_id, , drop = FALSE] * 1000
  # Anchor point chosen so the simulated domain falls entirely EAST of the
  # real example_geographic_filter()'s north-west exclusion box
  # (longitude < 34.80 & latitude > -0.46), which is calibrated to the real
  # Kisumu HDSS survey geography, not to this synthetic domain. Using the
  # real survey's own anchor (~34.75) here would cause that filter to drop
  # roughly two-thirds of the simulated cohort as a silent side effect.
  # lon0 = 34.95 keeps longitude > 34.80 for the domain width used below.
  lon0 <- 34.95; lat0 <- -0.45
  m_per_deg_lon <- 111000 * cos(lat0 * pi / 180)
  m_per_deg_lat <- 111000
  longitude <- lon0 + loc_m[, 1L] / m_per_deg_lon
  latitude  <- lat0 + loc_m[, 2L] / m_per_deg_lat

  data.frame(
    id        = seq_len(n_individual),
    age       = age,
    long_m    = loc_m[, 1L],
    lat_m     = loc_m[, 2L],
    longitude = longitude,
    latitude  = latitude,
    ama_norm  = ama_norm,
    msp_norm  = msp_norm
  )
}

# Minimal bivariate-normal sampler with unit variances and correlation rho,
# used only to avoid an MASS dependency for one call site above.
MASS_mvrnorm_fallback <- function(n, rho) {
  z1 <- stats::rnorm(n)
  z2 <- rho * z1 + sqrt(1 - rho^2) * stats::rnorm(n)
  cbind(z1, z2)
}

if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  sim_redhot <- simulate_redhot()
  dir.create("data", showWarnings = FALSE, recursive = TRUE)
  saveRDS(sim_redhot, file.path("data", "sim_redhot.rds"))
  cat(sprintf(
    "Simulated %d individuals at %d unique locations -> data/sim_redhot.rds\n",
    nrow(sim_redhot), length(unique(paste(sim_redhot$long_m, sim_redhot$lat_m)))
  ))
}
