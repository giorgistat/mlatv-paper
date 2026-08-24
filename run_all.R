# Runs the full pipeline (00-06) in order. Script 07 (simulation study) is a
# stub and is NOT run automatically -- see R/07_simulation_study.R and the
# README for why.
#
# Usage: Rscript run_all.R   (from the repository root)

repo_root <- getwd()
scripts <- c(
  "R/00_simulate_data.R",
  "R/01_univariate_ama1.R",
  "R/02_univariate_msp1.R",
  "R/03_bivariate_nonspatial.R",
  "R/04_bivariate_spatial_joint.R",
  "R/05_spatial_field_maps.R",
  "R/06_cross_validation_hdr.R"
)

for (s in scripts) {
  cat(sprintf("\n=========== Running %s ===========\n", s))
  t0 <- proc.time()
  source(file.path(repo_root, s), local = new.env(parent = globalenv()))
  cat(sprintf("--- %s finished in %.1fs ---\n", s, (proc.time() - t0)[["elapsed"]]))
}
cat("\nAll scripts finished. See output/rds/ and output/plot/.\n")
