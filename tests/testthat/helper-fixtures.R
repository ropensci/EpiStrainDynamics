# tests/testthat/helper-fixtures.R
# Downloads pre-fitted model fixtures from GitHub Release via piggyback.
# Files that use these: test-diagnose_model.R, test-metrics.R,
# test-plot.R, test-extended-regression.R

skip_if_no_extended_fixtures <- function() {
  skip_if_not(
    exists("extended_fixtures") && !is.null(extended_fixtures$fit_rw_single),
    "Extended fixtures not available - fixture download may have failed."
  )
}

if (!requireNamespace("piggyback", quietly = TRUE)) {
  skip("Package 'piggyback' is required to download test fixtures. Install with install.packages('piggyback').")
}

# Bridge GITHUB_TOKEN to GITHUB_PAT for piggyback authentication
if (nchar(Sys.getenv("GITHUB_TOKEN")) > 0 && nchar(Sys.getenv("GITHUB_PAT")) == 0) {
  Sys.setenv(GITHUB_PAT = Sys.getenv("GITHUB_TOKEN"))
}

message("EpiStrainDynamics: downloading test fixtures via piggyback...")

fixture_dir <- file.path(tempdir(), "testfixtures")
dir.create(fixture_dir, showWarnings = FALSE)

download_ok <- tryCatch(
  {
    piggyback::pb_download(
      repo      = "ropensci/epistraindynamics",
      tag       = "test-fixtures-v2",
      dest      = fixture_dir,
      overwrite = TRUE
    )
    TRUE
  },
  error = function(e) {
    message("EpiStrainDynamics: fixture download failed - ", conditionMessage(e))
    FALSE
  }
)

if (!download_ok) {
  skip("Fixture download from GitHub Release failed. Check network access and GitHub PAT.")
}

# Regular fixtures (small data subsets, minimal MCMC settings)
fit_rw_single <<- readRDS(file.path(fixture_dir, "fit_rw_single.rds"))
fit_ps_single <<- readRDS(file.path(fixture_dir, "fit_ps_single.rds"))
fit_rw_multi <<- readRDS(file.path(fixture_dir, "fit_rw_multi.rds"))
fit_ps_multi <<- readRDS(file.path(fixture_dir, "fit_ps_multi.rds"))
fit_rw_subtyped <<- readRDS(file.path(fixture_dir, "fit_rw_subtyped.rds"))
fit_ps_subtyped <<- readRDS(file.path(fixture_dir, "fit_ps_subtyped.rds"))
fit_rw_single_dow <<- readRDS(file.path(fixture_dir, "fit_rw_single_dow.rds"))
gi_simple <<- readRDS(file.path(fixture_dir, "gi_simple.rds"))
inc_single <<- readRDS(file.path(fixture_dir, "inc_single.rds"))
inc_multi <<- readRDS(file.path(fixture_dir, "inc_multi.rds"))
gr_single <<- readRDS(file.path(fixture_dir, "gr_single.rds"))
gr_multi <<- readRDS(file.path(fixture_dir, "gr_multi.rds"))
rt_single <<- readRDS(file.path(fixture_dir, "rt_single.rds"))
rt_multi <<- readRDS(file.path(fixture_dir, "rt_multi.rds"))
prop <<- readRDS(file.path(fixture_dir, "prop.rds"))

# Extended fixtures (full package data, thorough MCMC settings)
extended_fixtures <- new.env(parent = emptyenv())
extended_fixtures$fit_rw_single <- readRDS(file.path(fixture_dir, "ext_fit_rw_single.rds"))
extended_fixtures$fit_ps_single <- readRDS(file.path(fixture_dir, "ext_fit_ps_single.rds"))
extended_fixtures$fit_rw_multi <- readRDS(file.path(fixture_dir, "ext_fit_rw_multi.rds"))
extended_fixtures$fit_ps_multi <- readRDS(file.path(fixture_dir, "ext_fit_ps_multi.rds"))
extended_fixtures$fit_rw_subtyped <- readRDS(file.path(fixture_dir, "ext_fit_rw_subtyped.rds"))
extended_fixtures$fit_ps_subtyped <- readRDS(file.path(fixture_dir, "ext_fit_ps_subtyped.rds"))

message("EpiStrainDynamics: all fixtures ready.")
