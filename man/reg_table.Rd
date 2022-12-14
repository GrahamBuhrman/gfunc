% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_functions.R
\name{reg_table}
\alias{reg_table}
\title{Generate summary table for linear and general linear models.}
\usage{
reg_table(
  model,
  hide = NULL,
  adjust.method = "none",
  var.cov = NULL,
  exp.est = FALSE
)
}
\arguments{
\item{model}{Model object, of class "lm" or "glm".}

\item{hide}{Character vector, if not NULL, indicates which model statistics to be hidden in the table. Possible statistics include: "term", "estimate", "std.error", "conf.level", "conf.low", "conf.high", "statistic", "df.error", "p.value".}

\item{adjust.method}{= Character vector, if not NULL, indicates the method to adjust p-values. See stats::p.adjust() for details. Further possible adjustment methods are "tukey", "scheffe", "sidak" and "none" to explicitly disable adjustment for emmGrid objects (from emmeans).}

\item{var.cov}{= variance-covariance matrix used to compute uncertainty estimates (e.g., for robust standard errors). This argument accepts a covariance matrix, a function which returns a covariance matrix, or a string which identifies the function to be used to compute the covariance matrix.}

\item{exp.est}{= Logical, indicating whether or not to exponentiate the coefficients (and related confidence intervals). This is typical for logistic regression, or more generally speaking, for models with log or logit links. It is also recommended to use exponentiate = TRUE for models with log-transformed response values. Note: Delta-method standard errors are also computed (by multiplying the standard errors by the transformed coefficients). This is to mimic behaviour of other software packages, such as Stata, but these standard errors poorly estimate uncertainty for the transformed coefficient. The transformed confidence interval more clearly captures this uncertainty. For compare_parameters(), exponentiate = "nongaussian" will only exponentiate coefficients from non-Gaussian families.}
}
\value{
A gt table object of the model summary.
}
\description{
Given a multivariate linear or general linear regression model, create a table that summarizes the model.
}
