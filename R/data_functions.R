
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
  #' @param model Model object, of class "lm"
  #' @param line_color A color, quoted, to be used as the color of plotted lines
  #' @param line_size An integer, sets the size of plotted lines
  #' @param point_color A color, quoted, to be used as the color of plotted outliers
  #' @param point_shape An integer, sets the shape of plotted points
  #' @param point_alpha A floating number between 0 and 1, sets the opacity of plotted points
  #'
  #' @return A list of ggplot plot objects in the following order: Residual vs Fit, QQ Plot, Scale-Location, Leverage Plot, grid of all plots
  #' @importFrom rlang .data
  #' @export
  #'

  # get metrics of linear regression model

  model_metrics <-
    broom::augment(model) %>%
    dplyr::select(cooks_d = ".cooksd",
                  fitted = ".fitted",
                  hat = ".hat",
                  resid = ".resid",
                  std_resid = ".std.resid") %>%
    dplyr::mutate(cooks_cutoff = dplyr::case_when(cooks_d <= 1 ~ 0,
                                                  cooks_d > 1 ~ 1))

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
    ggplot2::ylim(min(resid),
                  max(resid)) +
    ggplot2::xlim(min(fitted),
                  max(fitted)) +
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
                                 y = sqrt(std_resid))) +
    ggplot2::geom_jitter(shape = point_shape,
                         alpha = point_alpha) +
    ggplot2::geom_smooth(se = FALSE,
                         color = line_color,
                         size = line_size) +
    ggplot2::ylim(min(sqrt(std_resid)),
                  max(sqrt(std_resid))) +
    ggplot2::xlim(min(fitted),
                  max(fitted)) +
    ggplot2::labs(x = "Fitted Value",
                  y = "Sqrt of Standardized Residual",
                  title = "Scale-Location") +
    ggplot2::theme_bw()

  # leverage plot

  leverage_plot <-
    ggplot2::ggplot(data = model_metrics,
                    ggplot2::aes(x = hat,
                                 y = std_resid,
                                 color = cooks_cutoff)) +
    ggplot2::geom_jitter(shape = point_shape,
                         alpha = point_alpha) +
    ggplot2::geom_smooth(se = FALSE,
                         color = line_color,
                         size = line_size) +
    ggplot2::ylim(min(std_resid),
                  max(std_resid)) +
    ggplot2::xlim(min(hat),
                  max(hat)) +
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

reg_table <- function(model, hide = NULL, adjust_method = "none", var_cov = NULL, exp_est = FALSE) {

  #' Generate MLR results table
  #'
  #' Given a multivariate linear regression model and the labels for parameters, create a results table
  #'
  #' @param model Model object, of class "lm" or "glm".
  #' @param hide Character vector, if not NULL, indicates which model statistics to be hidden in the table. Possible statistics include: "term", "estimate", "std.error", "conf.level", "conf.low", "conf.high", "statistic", "df.error", "p.value".
  #' @param adjust_method = Character vector, if not NULL, indicates the method to adjust p-values. See stats::p.adjust() for details. Further possible adjustment methods are "tukey", "scheffe", "sidak" and "none" to explicitly disable adjustment for emmGrid objects (from emmeans).
  #' @param var_cov = variance-covariance matrix used to compute uncertainty estimates (e.g., for robust standard errors). This argument accepts a covariance matrix, a function which returns a covariance matrix, or a string which identifies the function to be used to compute the covariance matrix.
  #' @param exp_est = Logical, indicating whether or not to exponentiate the coefficients (and related confidence intervals). This is typical for logistic regression, or more generally speaking, for models with log or logit links. It is also recommended to use exponentiate = TRUE for models with log-transformed response values. Note: Delta-method standard errors are also computed (by multiplying the standard errors by the transformed coefficients). This is to mimic behaviour of other software packages, such as Stata, but these standard errors poorly estimate uncertainty for the transformed coefficient. The transformed confidence interval more clearly captures this uncertainty. For compare_parameters(), exponentiate = "nongaussian" will only exponentiate coefficients from non-Gaussian families.
  #'
  #' @return A gt table object of the model summary.
  #' @importFrom rlang .data
  #' @export
  #'

  # generate results table for multivariate linear regression model

    if (class(model)[1] == "lm"){

      fit_stats <-
        broom::glance(model) %>%
        dplyr::select(`R<sup>2</sup>` = .data$r.squared,
                      `Adj R<sup>2</sup>` = .data$adj.r.squared,
                      .data$AIC,
                      .data$BIC) %>%
        dplyr::mutate(`f<sup>2</sup>` = .data$`R<sup>2</sup>` / (1 - .data$`R<sup>2</sup>`)) %>%
        dplyr::mutate_all(function(x) gtsummary::style_sigfig(x, digits = 3)) %>%
        {paste(names(.data), .data, sep = " = ", collapse = "; ")}

      table <-
        parameters::model_parameters(model = model,
                                     vcov = var_cov,
                                     p_adjust = adjust_method) %>%
        insight::standardize_names("broom") %>%
        insight::standardize_column_order("easystats") %>%
        tibble::tibble() %>%
        gt::gt() %>%
        gt::fmt_number(
          columns = c("estimate",
                      "std.error",
                      "conf.level",
                      "conf.low",
                      "conf.high",
                      "statistic",
                      "df.error",
                      "p.value"
                      ),
          decimals = 3
        ) %>%
        gt::tab_style(
          style = list(
            gt::cell_text(weight = "bold")
          ),
          locations = list(
            gt::cells_body(
              columns = "p.value",
              rows = "p.value" < 0.05
            ),
            gt::cells_column_labels()
          )
        ) %>%
        gt::tab_source_note(gt::html(fit_stats))

      if (!is.null(hide)) {

        table <-
          table %>%
          gt::cols_hide(columns = hide)

      }

    } else if (class(model)[1] == "glm") {

      foo <-
        data.frame(`Tjur's R<sup>2</sup>` = performance::r2_tjur(model)) %>%
        tibble::tibble()

      fit_stats <-
        foo %>%
        dplyr::mutate_all(function(x) gtsummary::style_sigfig(x, digits = 3)) %>%
        {paste(names(.data), .data, sep = " = ", collapse = "; ")}

      table <-
        parameters::model_parameters(model = model,
                                     exponentiate = exp_est,
                                     vcov = var_cov,
                                     p_adjust = adjust_method) %>%
        insight::standardize_names("broom") %>%
        insight::standardize_column_order("easystats") %>%
        tibble::tibble() %>%
        gt::gt() %>%
        gt::fmt_number(
          columns = c("estimate",
                      "std.error",
                      "conf.level",
                      "conf.low",
                      "conf.high",
                      "statistic",
                      "df.error",
                      "p.value"
                      ),
          decimals = 3
        ) %>%
        gt::tab_style(
          style = list(
            gt::cell_text(weight = "bold")
          ),
          locations = list(
            gt::cells_body(
              columns = "p.value",
              rows = "p.value" < 0.05
            ),
            gt::cells_column_labels()
          )
        ) %>%
        gt::tab_source_note(gt::html(fit_stats))

      if (!is.null(hide)) {

        table <-
          table %>%
          gt::cols_hide(columns = hide)

      }

    }

    return(table)

  }
