#' Mock data for cross-sectional studies
#'
#' Data generated for cross-sectional studies
#'
#' @name cross
#' @format `cross` and `small_cross`
#' A vector with 4 entries of matrix \eqn{y}, total count is 38000 for `cross` and 38 for `small_cross`:
#' \describe{
#'   \item{element 1}{\eqn{y_{1,1}}}
#'   \item{element 2}{\eqn{y_{1,2}}}
#'   \item{element 3}{\eqn{y_{2,1}}}
#'   \item{element 4}{\eqn{y_{2,2}}}
#' }
"cross"

#' @rdname cross
"small_cross"


#' Mock data for incomplete binary data
#'
#' Data generated for inference with incomplete binary data
#'
#' @name sat
#' @format `sat` and `small_sat`
#' A list of x, y, and r. The number of observations \eqn{n} is 3000 for `sat` and 100 for `small_sat`.
#' The number of predictors \eqn{p} in \eqn{X} is 4 for `sat` and 2 for `small_sat.`
#' \describe{
#'   \item{x}{A design matrix of X, binary. \eqn{2^p} columns and \eqn{n} rows.}
#'   \item{y}{A binary response variable. It contains missing entries.}
#'   \item{r}{A binary indicator variable. \eqn{R_i = 1} when \eqn{Y_i} is observed.}
#' }
"sat"

#' @rdname sat
"small_sat"


#' Mock data for incomplete count data
#'
#' Data generated for inference with incomplete count data
#'
#' @name count
#' @format `count` and `small_count`
#' A list of y and r. The number of observations \eqn{n} is 3000 for `count` and 100 for `small_count`.
#' \describe{
#'   \item{y}{A counr variable. It contains missing entries.}
#'   \item{r}{A binary indicator variable. \eqn{R_i = 1} when \eqn{Y_i} is observed.}
#' }
"count"

#' @rdname count
"small_count"

