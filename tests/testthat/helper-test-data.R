# Helper functions for creating test data and models using package datasets

# Create environment for storing mock objects
test_env <- new.env(parent = emptyenv())

# Helper function to check if package data exists
check_package_data <- function() {
  if (!exists("sarscov2") || !exists("influenza")) {
    skip("Package data (sarscov2, influenza) not available")
  }
}

# Helper function for testing method constructor structure
test_method_structure <- function(result, expected_method, expected_params = NULL) {
  # Test basic structure
  expect_type(result, "list")
  expect_s3_class(result, "EpiStrainDynamics.method")
  expect_equal(result$method, expected_method)
  expect_true("method" %in% names(result))
  expect_type(result$method, "character")

  # Test parameters if expected
  if (!is.null(expected_params)) {
    expect_true("model_params" %in% names(result))
    expect_type(result$model_params, "list")
    expect_named(result$model_params, names(expected_params))

    for (param_name in names(expected_params)) {
      expect_equal(result$model_params[[param_name]], expected_params[[param_name]])
    }
  } else {
    expect_length(result, 1)
  }
}

# Create standard test models using package datasets with various parameter combinations
create_test_models <- function() {
  check_package_data()

  # Base pathogen structures (reused across models)
  single_struct <- single(
    data = sarscov2,
    case_timeseries = "cases",
    time = "date"
  )

  multiple_struct <- multiple(
    data = sarscov2,
    case_timeseries = "cases",
    time = "date",
    component_pathogen_timeseries = c("alpha", "delta", "omicron", "other")
  )

  subtyped_struct <- subtyped(
    data = influenza,
    case_timeseries = "ili",
    time = "week",
    parent_timeseries = "inf_A",
    subtyped_timeseries = c("inf_H3N2", "inf_H1N1"),
    other_pathogen_timeseries = c("inf_B", "other")
  )

  list(
    # Basic models with default parameters
    rw_single = construct_model(
      pathogen_structure = single_struct,
      method = random_walk()
    ),
    ps_single = construct_model(
      pathogen_structure = single_struct,
      method = p_spline()
    ),
    rw_multiple = construct_model(
      pathogen_structure = multiple_struct,
      method = random_walk()
    ),
    ps_multiple = construct_model(
      pathogen_structure = multiple_struct,
      method = p_spline()
    ),
    rw_subtyped = construct_model(
      pathogen_structure = subtyped_struct,
      method = random_walk()
    ),
    ps_subtyped = construct_model(
      pathogen_structure = subtyped_struct,
      method = p_spline()
    ),

    # Models with custom smoothing structure (shared)
    rw_multiple_shared_smooth = construct_model(
      pathogen_structure = multiple_struct,
      method = random_walk(),
      smoothing_params = smoothing_structure("shared", tau_mean = 0, tau_sd = 1)
    ),
    ps_multiple_shared_smooth = construct_model(
      pathogen_structure = multiple_struct,
      method = p_spline(),
      smoothing_params = smoothing_structure("shared", tau_mean = 0, tau_sd = 1)
    ),

    # Models with independent smoothing
    rw_multiple_indep_smooth = construct_model(
      pathogen_structure = multiple_struct,
      method = random_walk(),
      smoothing_params = smoothing_structure("independent",
        tau_mean = c(0, 0.1, 0.3, 0),
        tau_sd = rep(1, 4)
      )
    ),
    ps_subtyped_indep_smooth = construct_model(
      pathogen_structure = subtyped_struct,
      method = p_spline(),
      smoothing_params = smoothing_structure("independent",
        tau_mean = c(0, 0, 0.2, 0.1),
        tau_sd = rep(0.5, 4)
      )
    ),

    # Models with correlated smoothing
    rw_multiple_corr_smooth = construct_model(
      pathogen_structure = multiple_struct,
      method = random_walk(),
      smoothing_params = smoothing_structure("correlated")
    ),
    ps_subtyped_corr_smooth = construct_model(
      pathogen_structure = subtyped_struct,
      method = p_spline(),
      smoothing_params = smoothing_structure("correlated")
    ),

    # Models with custom dispersion parameters
    rw_multiple_custom_disp = construct_model(
      pathogen_structure = multiple_struct,
      method = random_walk(),
      dispersion_params = dispersion_structure(phi_mean = 2.0, phi_sd = 0.5)
    ),
    ps_subtyped_custom_disp = construct_model(
      pathogen_structure = subtyped_struct,
      method = p_spline(),
      dispersion_params = dispersion_structure(phi_mean = 1.5, phi_sd = 0.75)
    ),

    # Models with pathogen noise enabled
    rw_multiple_noise = construct_model(
      pathogen_structure = multiple_struct,
      method = random_walk(),
      pathogen_noise = TRUE
    ),
    ps_subtyped_noise = construct_model(
      pathogen_structure = subtyped_struct,
      method = p_spline(),
      pathogen_noise = TRUE
    ),

    # Models with day-of-week effects
    rw_single_dow = construct_model(
      pathogen_structure = single_struct,
      method = random_walk(),
      dow_effect = TRUE
    ),
    ps_multiple_dow = construct_model(
      pathogen_structure = multiple_struct,
      method = p_spline(),
      dow_effect = TRUE
    ),

    # Combined features: custom smoothing + dispersion + noise + dow
    rw_multiple_full = construct_model(
      pathogen_structure = multiple_struct,
      method = random_walk(),
      smoothing_params = smoothing_structure("independent",
        tau_mean = rep(0, 4),
        tau_sd = rep(1, 4)
      ),
      dispersion_params = dispersion_structure(phi_mean = 2.0, phi_sd = 0.5),
      pathogen_noise = TRUE,
      dow_effect = TRUE
    ),
    ps_subtyped_full = construct_model(
      pathogen_structure = subtyped_struct,
      method = p_spline(),
      smoothing_params = smoothing_structure("correlated"),
      dispersion_params = dispersion_structure(phi_mean = 1.5, phi_sd = 0.3),
      pathogen_noise = TRUE,
      dow_effect = TRUE
    )
  )
}

# Helper function to get expected data lengths from package datasets
get_expected_data_lengths <- function() {
  check_package_data()

  list(
    sarscov2_length = length(sarscov2$cases),
    influenza_length = length(influenza$ili),
    sarscov2_pathogen_count = 4, # alpha, delta, omicron, other
    influenza_subtyped_count = 2, # H3N2, H1N1
    influenza_other_count = 2 # influenzaB, other
  )
}

# Helper function to get expected pathogen names
get_expected_pathogen_names <- function() {
  list(
    sarscov2_multiple = c("alpha", "delta", "omicron", "other"),
    influenza_subtyped = c("inf_H3N2", "inf_H1N1", "inf_B", "other")
  )
}
