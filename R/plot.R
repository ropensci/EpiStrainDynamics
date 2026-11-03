#' Plot metrics calculation outputs
#'
#' S3 methods for plotting metrics calculated from the output sof either
#' [incidence()], [growth_rate()], [Rt()], or [proportion()].
#'
#' @param x Metrics calculation output of class `EpiStrainDynamics.metric`
#' @param xlab Time label for x axis, defaults to "Time"
#' @param ... Additional arguments passed to plot
#' @importFrom viridis viridis
#' @importFrom stats setNames
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_point geom_hline
#'  theme_bw theme scale_colour_manual scale_y_continuous sec_axis xlab ylab
#'  ggplot_build element_blank labs
#' @importFrom rlang .data
#'
#' @return ggplot2 plot output
#' @name plot
#' @export
#'
#' @srrstats {BS6.1} All calculated epidemiological metrics can be plotted with
#'   a default plot function.
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {TS5.0} default plot methods implemented
#' @srrstats {TS5.1, TS5.2, TS5.3} Time on x axis with units printed
#'
#' @examplesIf interactive()
#' mod <- construct_model(
#'   pathogen_structure = single(
#'     case_timeseries = sarscov2$cases,
#'     time = sarscov2$date
#'   ),
#'   method = random_walk()
#' )
#'
#' fit <- fit_model(mod)
#' gr <- growth_rate(mod)
#' plot(gr)
plot.incidence <- function(x, xlab = "Time", ...) {
  validate_class_inherits(x, "EpiStrainDynamics.metric")

  measure_df <- x$measure

  total_pathogens <- "Total"
  all_levels <- unique(measure_df$pathogen)
  other_levels <- setdiff(all_levels, total_pathogens)

  colors <- c(
    setNames("black", total_pathogens),
    setNames(viridis::viridis(length(other_levels), end = 0.9), other_levels)
  )

  tsbl <- x$constructed_model$validated_tsbl
  time_col <- tsibble::index_var(tsbl)
  input_data <- tsbl[, c(time_col, "case_timeseries")]

  ggplot2::ggplot(measure_df) +
    ggplot2::geom_line(ggplot2::aes(
      x = .data$time,
      y = .data$y,
      color = .data$pathogen
    )) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        x = .data$time,
        y = .data$y,
        ymin = .data$lb_50,
        ymax = .data$ub_50,
        fill = .data$pathogen
      ),
      alpha = 0.2
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        x = .data$time,
        y = .data$y,
        ymin = .data$lb_95,
        ymax = .data$ub_95,
        fill = .data$pathogen
      ),
      alpha = 0.2
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::geom_point(
      data = input_data,
      ggplot2::aes(
        x = .data[[time_col]],
        y = .data$case_timeseries
      ),
      size = 0.8
    ) +
    ggplot2::geom_line(
      data = input_data,
      ggplot2::aes(
        x = .data[[time_col]],
        y = .data$case_timeseries
      ),
      linewidth = 0.2
    ) +
    ggplot2::scale_colour_manual(
      values = colors,
      aesthetics = c("colour", "fill")
    ) +
    ggplot2::ylab("Modelled cases") +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::xlab(xlab)
}

#' @rdname plot
#' @export
plot.growth_rate <- function(x, xlab = "Time", ...) {
  validate_class_inherits(x, "EpiStrainDynamics.metric")

  measure_df <- x$measure

  total_pathogens <- "Total"
  all_levels <- unique(measure_df$pathogen)
  other_levels <- setdiff(all_levels, total_pathogens)

  colors <- c(
    setNames("black", total_pathogens),
    setNames(viridis::viridis(length(other_levels), end = 0.9), other_levels)
  )

  p <- ggplot2::ggplot(measure_df) +
    ggplot2::geom_line(ggplot2::aes(
      x = .data$time,
      y = .data$y,
      color = .data$pathogen
    )) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        x = .data$time,
        y = .data$y,
        ymin = .data$lb_50,
        ymax = .data$ub_50,
        fill = .data$pathogen
      ),
      alpha = 0.2
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        x = .data$time,
        y = .data$y,
        ymin = .data$lb_95,
        ymax = .data$ub_95,
        fill = .data$pathogen
      ),
      alpha = 0.2
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::scale_colour_manual(
      values = colors,
      aesthetics = c("colour", "fill")
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::ylab("Growth rate") +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::xlab(xlab)

  # Get the y-axis range from the built plot
  y_range <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y.range

  # Calculate appropriate breaks
  n_breaks <- 5
  breaks <- pretty(y_range, n = n_breaks)

  # Create labels
  times <- log(2) / abs(breaks)
  signs <- ifelse(breaks >= 0, "", "-")
  labels <- ifelse(abs(breaks) < 1e-6, Inf,
    paste0(signs, round(times, 1))
  )

  # Add secondary axis with calculated breaks and labels
  p + ggplot2::scale_y_continuous(
    sec.axis = ggplot2::sec_axis(~.,
      breaks = breaks,
      labels = labels,
      name = "Doubling(+) / Halving(-) time"
    )
  )
}

#' @rdname plot
#' @export
plot.Rt <- function(x, xlab = "Time", ...) {
  validate_class_inherits(x, "EpiStrainDynamics.metric")

  measure_df <- x$measure

  total_pathogens <- "Total"
  all_levels <- unique(measure_df$pathogen)
  other_levels <- setdiff(all_levels, total_pathogens)

  colors <- c(
    setNames("black", total_pathogens),
    setNames(viridis::viridis(length(other_levels), end = 0.9), other_levels)
  )

  ggplot2::ggplot(measure_df) +
    ggplot2::geom_line(ggplot2::aes(
      x = .data$time,
      y = .data$y,
      color = .data$pathogen
    )) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        x = .data$time,
        y = .data$y,
        ymin = .data$lb_50,
        ymax = .data$ub_50,
        fill = .data$pathogen
      ),
      alpha = 0.2
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        x = .data$time,
        y = .data$y,
        ymin = .data$lb_95,
        ymax = .data$ub_95,
        fill = .data$pathogen
      ),
      alpha = 0.2
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::scale_colour_manual(
      values = colors,
      aesthetics = c("colour", "fill")
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::ylab("Effective reproduction number") +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::xlab(xlab)
}

#' @rdname plot
#' @export
plot.proportion <- function(x, xlab = "Time", ...) {
  validate_class_inherits(x, "EpiStrainDynamics.metric")

  measure_df <- x$measure

  combos <- unique(measure_df$pathogen)
  colors <- if (length(combos) == 1) {
    setNames("black", combos)
  } else {
    setNames(viridis::viridis(length(combos), end = 0.9), combos)
  }

  p <- ggplot2::ggplot(measure_df) +
    ggplot2::geom_line(ggplot2::aes(
      x = .data$time,
      y = .data$y,
      color = .data$pathogen
    )) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        x = .data$time,
        y = .data$y,
        ymin = .data$lb_50,
        ymax = .data$ub_50,
        fill = .data$pathogen
      ),
      alpha = 0.2
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        x = .data$time,
        y = .data$y,
        ymin = .data$lb_95,
        ymax = .data$ub_95,
        fill = .data$pathogen
      ),
      alpha = 0.2
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::scale_colour_manual(
      values = colors,
      aesthetics = c("colour", "fill")
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::ylab("Modelled proportion of cases") +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::xlab(xlab)

  # Surface a non-default denominator: it has no other representation in the
  # plot (unlike the numerator, which already appears in the pathogen legend).
  all_pathogens <- unique(x$constructed_model$pathogen_names)
  if (!is.null(x$denominator_combination) &&
    !setequal(x$denominator_combination, all_pathogens)) {
    p <- p + ggplot2::labs(
      subtitle = paste(
        "Denominator:",
        paste(x$denominator_combination, collapse = ", ")
      )
    )
  }

  p
}

#' Plot convergence diagnostics over time
#'
#' Plots the R-hat or effective sample size of the log-incidence trajectory
#' (`log_incidence`) against time, with one line per pathogen, so that
#' convergence problems can be located in time and by pathogen. The
#' threshold used by [diagnose_model()] is shown as a dashed line.
#'
#' R-hat is shown between 1 and 1.2, so values far above the threshold run
#' off the top of the plot. Effective sample size is shown on a log scale,
#' covering the range of the values and the threshold, so that values near
#' the threshold stay readable even when others are in the thousands.
#'
#' Only the trajectory is plotted. Other parameters, including the spline
#' coefficients of p-spline models, are listed in the `parameters` table of
#' the [diagnose_model()] output.
#'
#' @param x Output of [diagnose_model()], of class
#'   `EpiStrainDynamics.diagnostics`
#' @param metric The diagnostic to plot: `"rhat"` (default) or `"n_eff"`
#' @param xlab Time label for x axis, defaults to "Time"
#' @param ... Not used
#'
#' @return ggplot2 plot output
#' @export
#'
#' @srrstats {BS6.1} convergence diagnostics can be plotted with a default
#'   plot function
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {TS5.0} default plot methods implemented
#' @srrstats {TS5.1, TS5.2, TS5.3} Time on x axis with units printed
#'
#' @examplesIf interactive()
#' mod <- construct_model(
#'   pathogen_structure = single(
#'     data = sarscov2,
#'     case_timeseries = "cases",
#'     time = "date"
#'   ),
#'   method = random_walk()
#' )
#' fit <- fit_model(mod)
#' diag <- diagnose_model(fit)
#' plot(diag)
#' plot(diag, metric = "n_eff")
plot.EpiStrainDynamics.diagnostics <- function(x,
                                               metric = c("rhat", "n_eff"),
                                               xlab = "Time",
                                               ...) {
  validate_class_inherits(x, "EpiStrainDynamics.diagnostics")
  metric <- rlang::arg_match(metric)

  trajectory <- x$parameters[x$parameters$parameter_group == "log_incidence", ]

  pathogen_levels <- unique(trajectory$pathogen)
  colors <- setNames(
    viridis::viridis(length(pathogen_levels), end = 0.9),
    pathogen_levels
  )

  if (metric == "rhat") {
    threshold <- x$rhat_threshold
    ylim <- c(1, 1.2)
    ylab <- "R-hat"
    y_scale <- ggplot2::scale_y_continuous()
  } else {
    threshold <- x$eff_sample_threshold
    ylim <- range(c(trajectory$n_eff, threshold), finite = TRUE)
    ylab <- "Effective sample size"
    y_scale <- ggplot2::scale_y_log10()
  }

  ggplot2::ggplot(trajectory) +
    ggplot2::geom_line(ggplot2::aes(
      x = .data$time,
      y = .data[[metric]],
      color = .data$pathogen
    )) +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed") +
    y_scale +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::ylab(ylab) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::xlab(xlab)
}
