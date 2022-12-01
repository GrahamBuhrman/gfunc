
mlr_table <- function(model, var_names){

  #' Generate MLR results table
  #'
  #' Given a multivariate linear regression model and the labels for parameters, create a results table
  #'
  #' @param model The MLR model whose results should be used
  #' @param var_names A named list with entries of the format parameter = "Parameter Label"
  #'
  #' @return A gt table object
  #' @export
  #'
  #' @examples
  #'
  #' # create MLR model
  #'
  #' m1 <- lm(mpg ~ wt + am)
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

  require(tidyverse)
  require(gtsummary)
  require(parameters)
  require(broom)

  # get model fit statistics

  fit_stats <-
    broom::glance(model) %>%
    select(`R<sup>2</sup>` = r.squared,
           `Adj R<sup>2</sup>` = adj.r.squared,
           AIC, BIC) %>%
    mutate(`f<sup>2</sup>` = `R<sup>2</sup>` / (1 - `R<sup>2</sup>`)) %>%
    mutate_all(function(x) style_sigfig(x, digits = 3)) %>%
    {paste(names(.), ., sep = " = ", collapse = "; ")}

  # generate results table

  table <-
    tbl_regression(
      x = model,
      tidy_fun = function(x, conf.int = T, conf.level = 0.95)
        parameters::model_parameters(
          model = x,
          vcov = "vcovHC") %>%
        insight::standardize_names(style = "broom") %>%
        tibble(),
      label = var_names,
      intercept = T) %>%
    bold_labels() %>%
    add_vif() %>%
    add_q(method = "bonferroni") %>%
    bold_p(q = T) %>%
    modify_column_unhide(column = std.error) %>%
    as_gt() %>%
    gt::tab_source_note(gt::html(fit_stats))

  return(table)

}
