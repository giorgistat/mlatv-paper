# Purpose: Posterior mean/SD maps of the spatial fields S_1 (AMA1), S_2 (MSP1)
#   from the joint spatial fit, and the population-level predictive targets
#   of Section 4.3 (E[T_k | S], E[Y_k | S]) evaluated at each unique location.
# Adapted from: example22_full_spatial_z_realdata.R, restructured to build
#   directly on 04's joint fit rather than a separate example19/example20
#   pipeline (see README "What's missing", item 2) -- this is NOT a
#   byte-for-byte port of example22.R.
# Output: output/rds/05_spatial_field_maps.rds,
#         output/plot/05_spatial_z_fields_mean.pdf,
#         output/plot/05_spatial_z_fields_sd.pdf

library(mlatv)
suppressPackageStartupMessages({ library(Matrix) })
source(file.path(if (exists("resolve_repo_root")) resolve_repo_root() else getwd(),
                  "R", "utils_helpers.R"))
repo_root <- resolve_repo_root()
plot_dir <- file.path(repo_root, "output", "plot")
rds_dir  <- file.path(repo_root, "output", "rds")

prev_path <- file.path(rds_dir, "04_bivariate_spatial_joint.rds")
if (!file.exists(prev_path)) stop("Run 04_bivariate_spatial_joint.R first.")
prev <- readRDS(prev_path)
fit_joint  <- prev$fit_joint
loc_unique <- prev$loc_unique
mesh       <- prev$mesh

# ---- Posterior mean and SD of S_1 (AMA1) and S_2 (MSP1) at unique locations ---
n_mesh <- mesh$n
n_loc  <- nrow(loc_unique)
H_hat  <- fit_joint$H_hat
A_uniq <- fit_joint$bsz_data$A   # projection matrix, mesh -> unique locations
if (is.null(A_uniq)) {
  A_uniq <- fmesher::fm_basis(mesh, loc = loc_unique)
}

cat("Computing posterior SD via Cholesky solve...\n")
L_H <- Matrix::Cholesky(H_hat, perm = TRUE, LDL = FALSE, super = FALSE)

rhs_1 <- rbind(Matrix::t(A_uniq), Matrix::Matrix(0, n_mesh, n_loc, sparse = TRUE))
rhs_2 <- rbind(Matrix::Matrix(0, n_mesh, n_loc, sparse = TRUE), Matrix::t(A_uniq))

V_1 <- Matrix::solve(L_H, rhs_1)
V_2 <- Matrix::solve(L_H, rhs_2)
var_S_1 <- as.numeric(Matrix::colSums(rhs_1 * V_1))
var_S_2 <- as.numeric(Matrix::colSums(rhs_2 * V_2))
sd_S_1  <- sqrt(pmax(var_S_1, 0))
sd_S_2  <- sqrt(pmax(var_S_2, 0))

field_summary <- data.frame(
  x = loc_unique[, 1L], y = loc_unique[, 2L],
  mean_S_1 = fit_joint$s_hat[, "S_A"], mean_S_2 = fit_joint$s_hat[, "S_M"],
  sd_S_1   = sd_S_1,                   sd_S_2   = sd_S_2
)
cat("Field summary (first 6 rows):\n")
print(head(round(field_summary, 4)))

# ---- Plots ----------------------------------------------------------------------
plot_field <- function(x, y, vals, label, palette = "Blue-Red 3", zero_center = TRUE) {
  if (zero_center) {
    rng <- max(abs(vals), na.rm = TRUE) + 1e-6
    brks <- seq(-rng, rng, length.out = 31L)
  } else {
    brks <- seq(0, max(vals, na.rm = TRUE) * 1.001, length.out = 31L)
  }
  cols <- grDevices::hcl.colors(length(brks) - 1L, palette = palette,
                                 rev = !zero_center)
  ix <- as.integer(cut(vals, breaks = brks, include.lowest = TRUE))
  graphics::plot(x, y, pch = 19, col = cols[ix], cex = 0.7,
                 xlab = "x (m)", ylab = "y (m)", main = label)
}

fig1 <- file.path(plot_dir, "05_spatial_z_fields_mean.pdf")
grDevices::pdf(fig1, width = 9, height = 4.2)
opar <- graphics::par(mfrow = c(1L, 2L), mar = c(4, 4, 2, 1) + 0.1)
plot_field(field_summary$x, field_summary$y, field_summary$mean_S_1,
           expression(E[S[1] ~ "|" ~ Y]))
plot_field(field_summary$x, field_summary$y, field_summary$mean_S_2,
           expression(E[S[2] ~ "|" ~ Y]))
graphics::par(opar)
grDevices::dev.off()
cat("Saved:", fig1, "\n")

fig2 <- file.path(plot_dir, "05_spatial_z_fields_sd.pdf")
grDevices::pdf(fig2, width = 9, height = 4.2)
opar <- graphics::par(mfrow = c(1L, 2L), mar = c(4, 4, 2, 1) + 0.1)
plot_field(field_summary$x, field_summary$y, field_summary$sd_S_1,
           expression(SD(S[1] ~ "|" ~ Y)), palette = "YlOrRd", zero_center = FALSE)
plot_field(field_summary$x, field_summary$y, field_summary$sd_S_2,
           expression(SD(S[2] ~ "|" ~ Y)), palette = "YlOrRd", zero_center = FALSE)
graphics::par(opar)
grDevices::dev.off()
cat("Saved:", fig2, "\n")

# ---- Save -------------------------------------------------------------------
saveRDS(
  list(field_summary = field_summary),
  file.path(rds_dir, "05_spatial_field_maps.rds")
)
cat("Saved: output/rds/05_spatial_field_maps.rds\n")
