
# Diagnostic Plots for Linear Regression Models ----

diagnostic_plots <- function(model,
                             line_color = "steelblue",
                             line_size = 1,
                             point_color = "blue",
                             point_shape = 1,
                             point_alpha = 0.6) {

  #' Generate diagnostic plots for linear regression model
  #'
  #' Given a linear regression model and aesthetic parameters, generate plots to check assumptions of linear regression
  #'
  #' @param model The linear regression model to create diagnostic plots for
  #' @param line_color The color of lines in the plots. Default is "steelblue"
  #' @param line_size The size of lines in the plots. Default is 1
  #' @param point_color The color of outliers in the plots. Default is "blue"
  #' @param point_shape The shape of points in the plots. Default is 1
  #' @param point_alpha The alpha of points in the plots. Default is 0.6
  #'
  #' @return A list of ggplot plot objects in the following order: Residual vs Fit, QQ Plot, Scale-Location, Leverage Plot, grid of all plots
  #' @importFrom rlang .data
  #' @export
  #'
  #' @examples
  #'
  #' # create linear regression model
  #'
  #' data(mtcars)
  #'
  #' m1 <- lm(mpg ~ wt + am, data = mtcars)
  #'
  #' # generate diagnostic plots
  #'
  #' diagnostic_plots(model = m1)
  #'

  # get metrics of linear regression model

  resid <- fitted <- std_resid <- hat <- cooks_d <- cooks_cutoff <- sqrt_std_resid <- NULL

  model_metrics <-
    broom::augment(model) %>%
    dplyr::mutate(resid = .data$.resid,
                  fitted = .data$.fitted,
                  std_resid = .data$.std.resid,
                  hat = .data$.hat,
                  cooks_d = .data$.cooksd) %>%
    dplyr::mutate(cooks_cutoff = dplyr::case_when(.data$cooks_d <= 1 ~ 0,
                                                  .data$cooks_d > 1 ~ 1),
                  sqrt_std_resid = sqrt(abs(.data$std_resid)))

  # get range summaries of metrics

  metric_stats <-
    model_metrics %>%
    dplyr::summarize(resid_range = range(.data$resid),
                     fitted_range = range(.data$fitted),
                     stdresid_range = range(.data$std_resid),
                     hat_range = range(.data$hat),
                     sqrt_stdresid_range = range(.data$sqrt_std_resid))

  # residual vs fit plot

  resid_fit <-
    ggplot2::ggplot(data = model_metrics,
                    ggplot2::aes(x = fitted,
                                 y = resid)) +
    ggplot2::geom_jitter(shape = point_shape,
                         alpha = point_alpha) +
    ggplot2::geom_smooth(se = FALSE,
                         color = line_color,
                         size =  line_size) +
    ggplot2::ylim(metric_stats$resid_range[1],
                  metric_stats$resid_range[2]) +
    ggplot2::xlim(metric_stats$fitted_range[1],
                  metric_stats$fitted_range[2]) +
    ggplot2::labs(x = "Fitted Value",
                  y = "Residual",
                  title = "Residual vs Fit") +
    ggplot2::theme_bw()

  # quantile-quantile plot

  qq_plot <-
    ggplot2::ggplot(data = model_metrics,
                    ggplot2::aes(x =  stats::qqnorm(std_resid)[[1]],
                                 y = std_resid)) +
    ggplot2::geom_jitter() +
    ggplot2::geom_abline(size = line_size) +
    ggplot2::labs(x = "Theoretical Quantiles",
                  y = "Standardized Residual",
                  title = "Normal Q-Q") +
    ggplot2::theme_bw()

  # scale-location plot

  scale_location <-
    ggplot2::ggplot(data = model_metrics,
                    ggplot2::aes(x = fitted,
                                 y = sqrt_std_resid)) +
    ggplot2::geom_jitter(shape = point_shape,
                         alpha = point_alpha) +
    ggplot2::geom_smooth(se = FALSE,
                         color = line_color,
                         size = line_size) +
    ggplot2::ylim(metric_stats$sqrt_stdresid_range[1],
                  metric_stats$sqrt_stdresid_range[2]) +
    ggplot2::xlim(metric_stats$fitted_range[1],
                  metric_stats$fitted_range[2]) +
    ggplot2::labs(x = "Fitted Value",
                  y = "Sqrt of Standardized Residual",
                  title = "Scale-Location") +
    ggplot2::theme_bw()

  # leverage plot

  leverage_plot <-
    ggplot2::ggplot(data = model_metrics,
                    ggplot2::aes(x = hat,
                                 y = std_resid,
                        color = factor(cooks_cutoff))) +
    ggplot2::geom_jitter(shape = point_shape,
                         alpha = point_alpha) +
    ggplot2::geom_smooth(se = FALSE,
                         color = line_color,
                         size = line_size) +
    ggplot2::ylim(metric_stats$stdresid_range[1],
                  metric_stats$stdresid_range[2]) +
    ggplot2::xlim(metric_stats$hat_range[1],
                  metric_stats$hat_range[2]) +
    ggplot2::scale_color_manual(values = c("black", point_color)) +
    ggplot2::labs(x = "Leverage",
                  y = "Standardized Residual",
                  title = "Standardized Residual vs Leverage") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")

  # grid of all plots

  all_plots <-
    ggpubr::ggarrange(
      resid_fit, qq_plot, scale_location, leverage_plot,
      ncol = 2,
      nrow = 2,
      legend = "none"
    )

  # create list of all plots to return


  plot_list <-list("resid-fit" = resid_fit,
                   "qq" = qq_plot,
                   "scale-location" = scale_location,
                   "leverage" = leverage_plot,
                   "all-plots" = all_plots)

  return(plot_list)

}

# MLR Results Table Function ----

mlr_table <- function(model, var_names) {

  #' Generate MLR results table
  #'
  #' Given a multivariate linear regression model and the labels for parameters, create a results table
  #'
  #' @param model The MLR model whose results should be used
  #' @param var_names A named list with entries of the format parameter = "Parameter Label"
  #'
  #' @return A gt table object
  #' @importFrom rlang .data
  #' @export
  #'
  #' @examples
  #'
  #' # create MLR model
  #'
  #' data(mtcars)
  #'
  #' m1 <- lm(mpg ~ wt + am, data = mtcars)
  #'
  #' # make named list of parameter labels
  #'
  #' parameter_names <- list(wt = "Weight", am = "Transmission")
  #'
  #' # generate results table
  #'
  #' mlr_table(model = m1, var_names = parameter_names)
  #'

  # generate results table for multivariate linear regression model

  # get model fit statistics

  fit_stats <-
    broom::glance(model) %>%
    dplyr::select(`R<sup>2</sup>` = .data$r.squared,
                  `Adj R<sup>2</sup>` = .data$adj.r.squared,
                  .data$AIC,
                  .data$BIC) %>%
    dplyr::mutate(`f<sup>2</sup>` = .data$`R<sup>2</sup>` / (1 - .data$`R<sup>2</sup>`)) %>%
    dplyr::mutate_all(function(x) gtsummary::style_sigfig(x, digits = 3)) %>%
    {paste(names(.data), .data, sep = " = ", collapse = "; ")}

  # generate results table

  table <-
    gtsummary::tbl_regression(
      x = model,
      tidy_fun = function(x, conf.int = TRUE, conf.level = 0.95) {
        parameters::model_parameters(
          model = x,
          vcov = "vcovHC") %>%
        insight::standardize_names(style = "broom") %>%
        tibble::tibble()
        },
      label = var_names,
      intercept = TRUE) %>%
    gtsummary::bold_labels() %>%
    gtsummary::add_vif() %>%
    gtsummary::add_q(method = "bonferroni") %>%
    gtsummary::bold_p(q = TRUE) %>%
    gtsummary::modify_column_unhide(column = .data$std.error) %>%
    gtsummary::as_gt() %>%
    gt::tab_source_note(gt::html(fit_stats))

  return(table)

}
