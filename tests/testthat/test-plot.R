# Tests for plot() functions
# Focus: plot structure, layers, aesthetics

# ==============================================================================
# TESTS: INPUT VALIDATION
# ==============================================================================

test_that("plot() validates input class", {
  bad_growth_rate <- structure(list(measure = data.frame()), class = "growth_rate")
  expect_error(
    plot(bad_growth_rate),
    "Input must be of class"
  )
})


# ==============================================================================
# TESTS: PLOT RETURNS CORRECT OBJECT TYPE
# ==============================================================================

test_that("plot.incidence() returns ggplot object", {
  p <- plot(inc_single)
  expect_s3_class(p, "ggplot")
})

test_that("plot.growth_rate() returns ggplot object", {
  p <- plot(gr_single)
  expect_s3_class(p, "ggplot")
})

test_that("plot.Rt() returns ggplot object", {
  p <- plot(rt_single)
  expect_s3_class(p, "ggplot")
})

test_that("plot.proportion() returns ggplot object", {
  p <- plot(prop)
  expect_s3_class(p, "ggplot")
})

# ==============================================================================
# TESTS: PLOT STRUCTURE AND LAYERS
# ==============================================================================

test_that("plot.incidence() has expected layers", {
  p <- plot(inc_single)

  # Get layer information
  layers <- vapply(p$layers, function(x) class(x$geom)[1], FUN.VALUE = character(1))

  # Should have line, ribbons, and points
  expect_true("GeomLine" %in% layers)
  expect_true("GeomRibbon" %in% layers)
  expect_true("GeomPoint" %in% layers)

  # Check number of ribbon layers (50% and 95% CIs)
  n_ribbons <- sum(layers == "GeomRibbon")
  expect_equal(n_ribbons, 2)
})

test_that("plot.growth_rate() has expected layers", {
  p <- plot(gr_single)

  layers <- vapply(p$layers, function(x) class(x$geom)[1], FUN.VALUE = character(1))

  # Should have line, ribbons, and horizontal line
  expect_true("GeomLine" %in% layers)
  expect_true("GeomRibbon" %in% layers)
  expect_true("GeomHline" %in% layers)

  # Check for secondary axis (doubling/halving time)
  expect_true(!is.null(p$scales$get_scales("y")))
})

test_that("plot.Rt() has expected layers", {
  p <- plot(rt_single)

  layers <- vapply(p$layers, function(x) class(x$geom)[1], FUN.VALUE = character(1))

  # Should have line, ribbons, and horizontal line at 1
  expect_true("GeomLine" %in% layers)
  expect_true("GeomRibbon" %in% layers)
  expect_true("GeomHline" %in% layers)
})

test_that("plot.proportion() has expected layers", {
  p <- plot(prop)

  layers <- vapply(p$layers, function(x) class(x$geom)[1], FUN.VALUE = character(1))

  # Should have line, ribbons, and horizontal line
  expect_true("GeomLine" %in% layers)
  expect_true("GeomRibbon" %in% layers)
  expect_true("GeomHline" %in% layers)
})

# ==============================================================================
# TESTS: CORRECT AESTHETICS MAPPING
# ==============================================================================

test_that("plot.incidence() maps aesthetics correctly", {
  p <- plot(inc_single)

  # Check main line layer aesthetics
  line_layer <- p$layers[[1]]
  line_mapping <- line_layer$mapping

  expect_equal(rlang::as_label(line_mapping$x), "time")
  expect_equal(rlang::as_label(line_mapping$y), "y")
  expect_equal(rlang::as_label(line_mapping$colour), "pathogen")
})

test_that("plot.growth_rate() maps aesthetics correctly", {
  p <- plot(gr_single)

  line_layer <- p$layers[[1]]
  line_mapping <- line_layer$mapping

  expect_equal(rlang::as_label(line_mapping$x), "time")
  expect_equal(rlang::as_label(line_mapping$y), "y")
  expect_equal(rlang::as_label(line_mapping$colour), "pathogen")
})

test_that("plot.Rt() maps aesthetics correctly", {
  p <- plot(rt_single)

  line_layer <- p$layers[[1]]
  line_mapping <- line_layer$mapping

  expect_equal(rlang::as_label(line_mapping$x), "time")
  expect_equal(rlang::as_label(line_mapping$y), "y")
  expect_equal(rlang::as_label(line_mapping$colour), "pathogen")
})

test_that("plot.proportion() maps aesthetics correctly", {
  p <- plot(prop)

  line_layer <- p$layers[[1]]
  line_mapping <- line_layer$mapping

  expect_equal(rlang::as_label(line_mapping$x), "time")
  expect_equal(rlang::as_label(line_mapping$y), "y")
  expect_equal(rlang::as_label(line_mapping$colour), "pathogen")
})

# ==============================================================================
# TESTS: MULTIPLE PATHOGENS COLORING
# ==============================================================================

test_that("plot.incidence() handles multiple pathogens with correct colors", {
  p <- plot(inc_multi)

  # Extract color scale
  color_scale <- p$scales$get_scales("colour")

  # Should have manual color scale
  expect_s3_class(color_scale, "ScaleDiscrete")

  # Check that "Total" is present and colored black
  expect_true("Total" %in% names(color_scale$palette(5)))
})

test_that("plot.growth_rate() handles multiple pathogens", {
  p <- plot(gr_multi)

  # Should have data for all pathogens
  pathogen_levels <- unique(p$data$pathogen)
  expect_true("Total" %in% pathogen_levels)
  expect_equal(length(pathogen_levels), 5) # 4 pathogens + Total
})

test_that("plot.Rt() handles multiple pathogens", {
  p <- plot(rt_multi)

  pathogen_levels <- unique(p$data$pathogen)
  expect_true("Total" %in% pathogen_levels)
  expect_equal(length(pathogen_levels), 5)
})

test_that("plot.proportion() uses distinct colors for combinations", {
  p <- plot(prop)

  # Check color scale
  color_scale <- p$scales$get_scales("colour")
  expect_s3_class(color_scale, "ScaleDiscrete")

  # Number of colors should match number of pathogen combinations
  n_combos <- length(unique(p$data$pathogen))
  expect_equal(n_combos, 4) # Default: each pathogen individually
})

# ==============================================================================
# TESTS: AXIS LABELS AND THEMES
# ==============================================================================

test_that("plot.incidence() has appropriate axis labels", {
  p <- plot(inc_single)

  # Check y-axis label
  expect_match(p$labels$y, "cases", ignore.case = TRUE)

  # Check theme
  expect_s3_class(p$theme, "theme")
})

test_that("plot.growth_rate() has appropriate axis labels", {
  p <- plot(gr_single)

  # Check y-axis label
  expect_match(p$labels$y, "Growth rate", ignore.case = TRUE)
})

test_that("plot.Rt() has appropriate axis labels", {
  p <- plot(rt_single)

  # Check y-axis label
  expect_match(p$labels$y, "reproduction", ignore.case = TRUE)
})

test_that("plot.proportion() has appropriate axis labels", {
  p <- plot(prop)

  # Check y-axis label
  expect_match(p$labels$y, "proportion", ignore.case = TRUE)
})

# ==============================================================================
# TESTS: DENOMINATOR SUBTITLE
# ==============================================================================

test_that("plot.proportion() has no subtitle for the default denominator", {
  skip_if_not(exists("fit_rw_multi"), "Cached fitted models not available")

  default_prop <- proportion(fit_rw_multi)
  p <- plot(default_prop)

  expect_null(p$labels$subtitle)
})

test_that("plot.proportion() shows a subtitle for a custom denominator", {
  skip_if_not(exists("fit_rw_multi"), "Cached fitted models not available")

  custom_prop <- proportion(fit_rw_multi,
    numerator_combination = "alpha",
    denominator_combination = c("alpha", "delta")
  )
  p <- plot(custom_prop)

  expect_match(p$labels$subtitle, "alpha")
  expect_match(p$labels$subtitle, "delta")
})

# ==============================================================================
# TESTS: REFERENCE LINES
# ==============================================================================

test_that("plot.growth_rate() includes horizontal line", {
  p <- plot(gr_single)

  # Find hline layer
  hline_layers <- which(vapply(p$layers, function(x) inherits(x$geom, "GeomHline"),
    FUN.VALUE = logical(1)
  ))
  expect_gt(length(hline_layers), 0)
})

test_that("plot.Rt() includes horizontal line", {
  p <- plot(rt_single)

  # Find hline layer
  hline_layers <- which(vapply(p$layers, function(x) inherits(x$geom, "GeomHline"),
    FUN.VALUE = logical(1)
  ))
  expect_gt(length(hline_layers), 0)
})

test_that("plot.proportion() includes horizontal line", {
  p <- plot(prop)

  # Find hline layer
  hline_layers <- which(vapply(p$layers, function(x) inherits(x$geom, "GeomHline"),
    FUN.VALUE = logical(1)
  ))
  expect_gt(length(hline_layers), 0)
})

# ==============================================================================
# TESTS: PLOT CAN BE RENDERED WITHOUT ERROR
# ==============================================================================

test_that("plot.incidence() can be rendered", {
  p <- plot(inc_single)

  # Should not error when building/printing
  expect_no_error(ggplot_build(p))

  # Suppress plot output during testing
  pdf(NULL)
  expect_no_error(print(p))
  dev.off()
})

test_that("plot.growth_rate() can be rendered", {
  p <- plot(gr_single)

  expect_no_error(ggplot_build(p))

  # Suppress plot output during testing
  pdf(NULL)
  expect_no_error(print(p))
  dev.off()
})

test_that("plot.Rt() can be rendered", {
  p <- plot(rt_single)

  expect_no_error(ggplot_build(p))

  # Suppress plot output during testing
  pdf(NULL)
  expect_no_error(print(p))
  dev.off()
})

test_that("plot.proportion() can be rendered", {
  p <- plot(prop)

  expect_no_error(ggplot_build(p))

  # Suppress plot output during testing
  pdf(NULL)
  expect_no_error(print(p))
  dev.off()
})

# ==============================================================================
# TESTS: CREDIBLE INTERVALS
# ==============================================================================

test_that("plot.incidence() includes both 50% and 95% credible intervals", {
  p <- plot(inc_single)

  # Find ribbon layers
  ribbon_layers <- which(vapply(p$layers, function(x) inherits(x$geom, "GeomRibbon"),
    FUN.VALUE = logical(1)
  ))
  expect_equal(length(ribbon_layers), 2)

  # Check that ribbons have ymin and ymax mappings
  for (i in ribbon_layers) {
    ribbon_mapping <- p$layers[[i]]$mapping
    expect_true(!is.null(ribbon_mapping$ymin))
    expect_true(!is.null(ribbon_mapping$ymax))
  }
})

test_that("plot.growth_rate() includes credible intervals", {
  p <- plot(gr_single)

  ribbon_layers <- which(vapply(p$layers, function(x) inherits(x$geom, "GeomRibbon"),
    FUN.VALUE = logical(1)
  ))
  expect_equal(length(ribbon_layers), 2)
})

test_that("plot.Rt() includes credible intervals", {
  p <- plot(rt_single)

  ribbon_layers <- which(vapply(p$layers, function(x) inherits(x$geom, "GeomRibbon"),
    FUN.VALUE = logical(1)
  ))
  expect_equal(length(ribbon_layers), 2)
})

test_that("plot.proportion() includes credible intervals", {
  p <- plot(prop)

  ribbon_layers <- which(vapply(p$layers, function(x) inherits(x$geom, "GeomRibbon"),
    FUN.VALUE = logical(1)
  ))
  expect_equal(length(ribbon_layers), 2)
})

# ==============================================================================
# TESTS: DATA POINTS (for incidence only)
# ==============================================================================

test_that("plot.incidence() includes original data points", {
  p <- plot(inc_single)

  # Should have point and line layers for original data
  layers <- vapply(p$layers, function(x) class(x$geom)[1],
    FUN.VALUE = character(1)
  )

  # Count points - should have at least one for raw data
  n_points <- sum(layers == "GeomPoint")
  expect_gte(n_points, 1)

  # Check that point layer uses different data (input_data)
  point_layer_idx <- which(layers == "GeomPoint")[1]
  point_data <- p$layers[[point_layer_idx]]$data
  expect_true("case_timeseries" %in% names(point_data))
})

# ==============================================================================
# TESTS: CONVERGENCE DIAGNOSTICS
# ==============================================================================

test_that("plot.EpiStrainDynamics.diagnostics() plots the trajectory against time", {
  invisible(capture.output(diag <- diagnose_model(fit_rw_subtyped)))
  trajectory <- diag$parameters[diag$parameters$parameter_group == "log_incidence", ]

  p <- plot(diag)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), nrow(trajectory))

  layers <- vapply(p$layers, function(x) class(x$geom)[1], FUN.VALUE = character(1))
  expect_true("GeomLine" %in% layers)
  expect_true("GeomHline" %in% layers)

  built <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(built$data[[1]]$colour)), 4)
  expect_equal(built$data[[2]]$yintercept, diag$rhat_threshold)
})

test_that("plot.EpiStrainDynamics.diagnostics() shows n_eff on a log scale covering the threshold", {
  invisible(capture.output(
    diag <- diagnose_model(fit_rw_subtyped, eff_sample_threshold = 500)
  ))

  p <- plot(diag, metric = "n_eff")
  expect_gte(max(p$coordinates$limits$y), diag$eff_sample_threshold)

  built <- ggplot2::ggplot_build(p)
  expect_equal(sort(built$data[[1]]$y), sort(log10(p$data$n_eff)))
  expect_equal(built$data[[2]]$yintercept, log10(500))
})

test_that("plot.EpiStrainDynamics.diagnostics() validates its arguments", {
  invisible(capture.output(diag <- diagnose_model(fit_rw_subtyped)))

  expect_error(plot(diag, metric = "bad"), "must be one of")
  expect_error(
    plot.EpiStrainDynamics.diagnostics(list()),
    "Input must be of class"
  )
})

# ==============================================================================
# OPTIONAL: VISUAL REGRESSION TESTS (requires vdiffr package)
# ==============================================================================

# Uncomment to use vdiffr for visual regression testing
# test_that("plot.incidence() visual appearance", {
#   skip_if_not_installed("vdiffr")
#   p <- plot(inc_single)
#   vdiffr::expect_doppelganger("incidence_single", p)
# })
#
# test_that("plot.growth_rate() visual appearance", {
#   skip_if_not_installed("vdiffr")
#   p <- plot(gr_single)
#   vdiffr::expect_doppelganger("growth_rate_single", p)
# })
#
# test_that("plot.Rt() visual appearance", {
#   skip_if_not_installed("vdiffr")
#   p <- plot(rt_single)
#   vdiffr::expect_doppelganger("rt_single", p)
# })
#
# test_that("plot.proportion() visual appearance", {
#   skip_if_not_installed("vdiffr")
#   p <- plot(prop)
#   vdiffr::expect_doppelganger("proportion_multi", p)
# })

# ==============================================================================
# TESTS: EDGE CASES
# ==============================================================================

test_that("plots handle missing pathogen names", {
  # Create modified metric without proper pathogen names
  bad_metric <- inc_single
  bad_metric$measure$pathogen <- NA

  # Should still create plot (though may look odd)
  expect_no_error(p <- plot(bad_metric))
  expect_s3_class(p, "ggplot")
})
