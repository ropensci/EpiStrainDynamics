# EpiStrainDynamics — test suite

## Running tests

Standard tests run with the usual `devtools::test()` or `R CMD check`. No
special configuration is needed. All tests, including regression tests against
pre-fitted models, run by default.

## How fixtures work

Rather than fitting Stan models during the test suite (which is slow), most
tests use fitted model fixtures that are pre-generated locally and stored in a
GitHub Release. They are downloaded automatically at the start of each test
session via [piggyback](https://docs.ropensci.org/piggyback/). The exceptions
are `test-fit_model.R` and `test-verbosity.R`, which test `fit_model()` itself
and so fit small models directly.

Two sets of fixtures are used, aided by `helper-fixtures.R`:

- **Regular fixtures**: small data subsets, minimal MCMC settings (500
  iterations, 1 chain). Used by `test-diagnose_model.R`, `test-metrics.R`,
  and `test-plot.R`.
- **Extended fixtures**: full package datasets, thorough MCMC settings (2000
  iterations, 4 chains). Used by `test-extended-regression.R`.

The helper file is sourced automatically by testthat before any tests run.
Fixtures live in memory for the duration of the test session and are discarded
afterwards.

## Regenerating fixtures

Fixtures only need to be regenerated when model structure changes in a way that
affects the fitted object format. To regenerate:
```r
devtools::load_all(helpers = FALSE)
source("inst/scripts/generate-test-fixtures.R")
```

`helpers = FALSE` stops `load_all()` running `helper-fixtures.R`, which would
otherwise try to download from the new tag before it has any files and fail.

This fits all regular and extended fixtures locally, saves them to
`inst/testfixtures/`, and uploads them to the GitHub Release. Bump the tag
(e.g. `test-fixtures-v2`) in the script and update the tag name in the helper 
file when doing so.

The upload needs the release to exist already. Create it first, marked as not
the latest release so it doesn't replace the package release on the repo page:
```sh
gh release create test-fixtures-v2 -R ropensci/EpiStrainDynamics \
  --latest=false --title test-fixtures-v2 --notes "Test fixtures"
```

## Adding new tests

Add test cases to the appropriate `test-*.R` file. If the test needs a new
fitted model, add it to `inst/scripts/generate-test-fixtures.R`, regenerate,
and load it in `helper-fixtures.R`.

## Requirements

| Requirement | Detail |
|---|---|
| **Platform** | No platform-specific requirements. Tests run on Linux, macOS, and Windows. |
| **Memory** | The fixture release is roughly 800 MB, and all fixtures are loaded into memory at once, so allow at least that much free RAM. |
| **Expected runtime** | A few minutes. Most tests use the fixtures; `test-fit_model.R` and `test-verbosity.R` fit small models directly. |
| **Artefacts** | No files are written to disk during tests. Fixtures are downloaded to a temporary directory and discarded afterwards. |
| **Dependencies** | The Stan models are fitted with rstan, which is installed as a package dependency; no separate Stan installation is needed. Downloading the fixtures needs `piggyback` and, to avoid GitHub rate limits, a `GITHUB_PAT` or `GITHUB_TOKEN`. |
