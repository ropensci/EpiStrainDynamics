# inst/scripts/generate-test-fixtures.R
#
# Generates and saves all test fixtures (regular and extended) to
# inst/testfixtures/ and uploads them to the GitHub Release via piggyback.
#
# Run this script locally whenever fixtures need to be regenerated (e.g. after
# model structure changes). Do NOT run on CI.
#
# Usage:
#   source("inst/scripts/generate-test-fixtures.R")
#
# Requirements:
#   - Run from the package root directory
#   - Package must be loaded: devtools::load_all(helpers = FALSE)
#     (helpers = FALSE skips helper-fixtures.R, which would try to download
#     from the new, still-empty release tag and fail)
#   - GitHub PAT with repo scope set in environment

library(piggyback)

# ------------------------------------------------------------------------------
# Checks
# ------------------------------------------------------------------------------

if (!dir.exists("inst/testfixtures")) {
  stop("inst/testfixtures does not exist - check your working directory")
}

if (!exists("sarscov2") || !exists("influenza")) {
  stop("Package data (sarscov2, influenza) not available - run devtools::load_all() first")
}

# Multi-chain fits run in parallel worker processes, which load the compiled
# models from the *installed* package rather than the one from load_all(). If
# the installed copy is out of date, those fixtures silently get the old
# models. Stop before saving if a fit lacks the current parameter names.
# Update `expected` if the Stan parameter names change.
check_fixture_pars <- function(fit, expected = "log_incidence") {
  missing <- setdiff(expected, fit$fit@model_pars)
  if (length(missing) > 0) {
    stop(
      deparse(substitute(fit)), " is missing parameter(s) ",
      paste(missing, collapse = ", "), ". The installed package is probably ",
      "out of date: run devtools::install() and rerun this script."
    )
  }
  invisible(fit)
}

# ------------------------------------------------------------------------------
# Regular fixtures
# (small data subsets, minimal MCMC settings for fast CI tests)
# ------------------------------------------------------------------------------

message("=== Generating regular fixtures ===")

sarscov2_subset  <- sarscov2[1:50, ]
influenza_subset <- influenza[1:40, ]

message("  Fitting: Random walk, single pathogen...")
fit_rw_single <- suppressWarnings(fit_model(
  construct_model(
    pathogen_structure = single(
      data = sarscov2_subset, case_timeseries = "cases", time = "date"
    ),
    method = random_walk()
  ),
  n_iter = 500, n_chain = 1, verbose = FALSE, seed = 123
))
check_fixture_pars(fit_rw_single)
saveRDS(fit_rw_single, "inst/testfixtures/fit_rw_single.rds")

message("  Fitting: P-spline, single pathogen...")
fit_ps_single <- suppressWarnings(fit_model(
  construct_model(
    pathogen_structure = single(
      data = sarscov2_subset, case_timeseries = "cases", time = "date"
    ),
    method = p_spline()
  ),
  n_iter = 500, n_chain = 1, verbose = FALSE, seed = 123
))
check_fixture_pars(fit_ps_single)
saveRDS(fit_ps_single, "inst/testfixtures/fit_ps_single.rds")

message("  Fitting: Random walk, multiple pathogens...")
fit_rw_multi <- suppressWarnings(fit_model(
  construct_model(
    pathogen_structure = multiple(
      data = sarscov2_subset, case_timeseries = "cases", time = "date",
      component_pathogen_timeseries = c("alpha", "delta", "omicron", "other")
    ),
    method = random_walk()
  ),
  n_iter = 500, n_chain = 1, verbose = FALSE, seed = 123
))
check_fixture_pars(fit_rw_multi)
saveRDS(fit_rw_multi, "inst/testfixtures/fit_rw_multi.rds")

message("  Fitting: P-spline, multiple pathogens...")
fit_ps_multi <- suppressWarnings(fit_model(
  construct_model(
    pathogen_structure = multiple(
      data = sarscov2_subset, case_timeseries = "cases", time = "date",
      component_pathogen_timeseries = c("alpha", "delta", "omicron", "other")
    ),
    method = p_spline()
  ),
  n_iter = 500, n_chain = 1, verbose = FALSE, seed = 123
))
check_fixture_pars(fit_ps_multi)
saveRDS(fit_ps_multi, "inst/testfixtures/fit_ps_multi.rds")

message("  Fitting: Random walk, subtyped pathogens...")
fit_rw_subtyped <- suppressWarnings(fit_model(
  construct_model(
    pathogen_structure = subtyped(
      data = influenza_subset, case_timeseries = "ili", time = "week",
      parent_timeseries = "inf_A",
      subtyped_timeseries = c("inf_H3N2", "inf_H1N1"),
      other_pathogen_timeseries = c("inf_B", "other")
    ),
    method = random_walk()
  ),
  n_iter = 500, n_chain = 1, verbose = FALSE, seed = 123
))
check_fixture_pars(fit_rw_subtyped)
saveRDS(fit_rw_subtyped, "inst/testfixtures/fit_rw_subtyped.rds")

message("  Fitting: P-spline, subtyped pathogens...")
fit_ps_subtyped <- suppressWarnings(fit_model(
  construct_model(
    pathogen_structure = subtyped(
      data = influenza_subset, case_timeseries = "ili", time = "week",
      parent_timeseries = "inf_A",
      subtyped_timeseries = c("inf_H3N2", "inf_H1N1"),
      other_pathogen_timeseries = c("inf_B", "other")
    ),
    method = p_spline()
  ),
  n_iter = 500, n_chain = 1, verbose = FALSE, seed = 123
))
check_fixture_pars(fit_ps_subtyped)
saveRDS(fit_ps_subtyped, "inst/testfixtures/fit_ps_subtyped.rds")

message("  Fitting: Random walk, single pathogen, day-of-week...")
fit_rw_single_dow <- suppressWarnings(fit_model(
  construct_model(
    pathogen_structure = single(
      data = sarscov2_subset, case_timeseries = "cases", time = "date"
    ),
    method = random_walk(),
    dow_effect = TRUE
  ),
  n_iter = 500, n_chain = 1, verbose = FALSE, seed = 123
))
check_fixture_pars(fit_rw_single_dow)
saveRDS(fit_rw_single_dow, "inst/testfixtures/fit_rw_single_dow.rds")

message("  Computing metrics...")
gi_simple <- function(x) ifelse(x == 0, 0, 4 * x * exp(-2 * x))
saveRDS(gi_simple,                          "inst/testfixtures/gi_simple.rds")
saveRDS(incidence(fit_rw_single, dow = FALSE), "inst/testfixtures/inc_single.rds")
saveRDS(incidence(fit_rw_multi,  dow = FALSE), "inst/testfixtures/inc_multi.rds")
saveRDS(growth_rate(fit_rw_single),            "inst/testfixtures/gr_single.rds")
saveRDS(growth_rate(fit_rw_multi),             "inst/testfixtures/gr_multi.rds")
saveRDS(Rt(fit_rw_single, gi_dist = gi_simple, tau_max = 7), "inst/testfixtures/rt_single.rds")
saveRDS(Rt(fit_rw_multi, gi_dist = gi_simple, tau_max = 7), "inst/testfixtures/rt_multi.rds")
saveRDS(proportion(fit_rw_multi),              "inst/testfixtures/prop.rds")

message("=== Regular fixtures done ===")

# ------------------------------------------------------------------------------
# Extended fixtures
# (full data, thorough MCMC settings for validation tests)
# ------------------------------------------------------------------------------

message("=== Generating extended fixtures ===")

message("  Fitting: Random walk, single pathogen...")
fit_rw_single_ext <- fit_model(
  construct_model(
    pathogen_structure = single(
      data = sarscov2, case_timeseries = "cases", time = "date"
    ),
    method = random_walk()
  ),
  n_iter = 2000, n_chain = 4, verbose = FALSE
)
check_fixture_pars(fit_rw_single_ext)
saveRDS(fit_rw_single_ext, "inst/testfixtures/ext_fit_rw_single.rds")

message("  Fitting: P-spline, single pathogen...")
fit_ps_single_ext <- fit_model(
  construct_model(
    pathogen_structure = single(
      data = sarscov2, case_timeseries = "cases", time = "date"
    ),
    method = p_spline()
  ),
  n_iter = 2000, n_chain = 4, verbose = FALSE
)
check_fixture_pars(fit_ps_single_ext)
saveRDS(fit_ps_single_ext, "inst/testfixtures/ext_fit_ps_single.rds")

message("  Fitting: Random walk, multiple pathogens...")
fit_rw_multi_ext <- fit_model(
  construct_model(
    pathogen_structure = multiple(
      data = sarscov2, case_timeseries = "cases", time = "date",
      component_pathogen_timeseries = c("alpha", "delta", "omicron", "other")
    ),
    method = random_walk()
  ),
  n_iter = 2000, n_chain = 4, verbose = FALSE
)
check_fixture_pars(fit_rw_multi_ext)
saveRDS(fit_rw_multi_ext, "inst/testfixtures/ext_fit_rw_multi.rds")

message("  Fitting: P-spline, multiple pathogens...")
fit_ps_multi_ext <- fit_model(
  construct_model(
    pathogen_structure = multiple(
      data = sarscov2, case_timeseries = "cases", time = "date",
      component_pathogen_timeseries = c("alpha", "delta", "omicron", "other")
    ),
    method = p_spline()
  ),
  n_iter = 2000, n_chain = 4, verbose = FALSE
)
check_fixture_pars(fit_ps_multi_ext)
saveRDS(fit_ps_multi_ext, "inst/testfixtures/ext_fit_ps_multi.rds")

message("  Fitting: Random walk, subtyped pathogens...")
fit_rw_subtyped_ext <- fit_model(
  construct_model(
    pathogen_structure = subtyped(
      data = influenza, case_timeseries = "ili", time = "week",
      parent_timeseries = "inf_A",
      subtyped_timeseries = c("inf_H3N2", "inf_H1N1"),
      other_pathogen_timeseries = c("inf_B", "other")
    ),
    method = random_walk()
  ),
  n_iter = 2000, n_chain = 4, verbose = FALSE
)
check_fixture_pars(fit_rw_subtyped_ext)
saveRDS(fit_rw_subtyped_ext, "inst/testfixtures/ext_fit_rw_subtyped.rds")

message("  Fitting: P-spline, subtyped pathogens...")
fit_ps_subtyped_ext <- fit_model(
  construct_model(
    pathogen_structure = subtyped(
      data = influenza, case_timeseries = "ili", time = "week",
      parent_timeseries = "inf_A",
      subtyped_timeseries = c("inf_H3N2", "inf_H1N1"),
      other_pathogen_timeseries = c("inf_B", "other")
    ),
    method = p_spline()
  ),
  n_iter = 2000, n_chain = 4, verbose = FALSE
)
check_fixture_pars(fit_ps_subtyped_ext)
saveRDS(fit_ps_subtyped_ext, "inst/testfixtures/ext_fit_ps_subtyped.rds")

message("=== Extended fixtures done ===")

# ------------------------------------------------------------------------------
# Upload all fixtures to GitHub Release via piggyback
# ------------------------------------------------------------------------------

message("=== Uploading fixtures to GitHub Release ===")

pb_upload(
  file      = list.files("inst/testfixtures", full.names = TRUE),
  repo      = "ropensci/epistraindynamics",
  tag       = "test-fixtures-v2"
)

message("=== All done. Fixtures uploaded to test-fixtures-v2 ===")
