#' Sparse Shooting S-Estimator
#'
#' Computes the sparse shooting S-estimator for robust regression in
#' high-dimensional settings with cellwise contamination. A BIC criterion
#' based on the residual M-scale is used to select the sparsity parameter
#' from a data-driven grid.
#'
#' @param x An `n x p` numeric matrix of predictors.
#' @param y A numeric vector of length `n` (the response).
#' @param k Tuning constant for Tukey's biweight rho function. Default is
#'   `3.420`, which yields 85% efficiency under normality.
#' @param maxIteration Maximum number of outer shooting loop iterations.
#'   Default is `100`.
#' @param tol Convergence tolerance for the outer shooting loop. Default is
#'   `1e-2`.
#' @param betaEst Optional `p x 1` matrix of initial regression coefficient
#'   estimates. If `NULL` (default), initialised internally via an MM-estimator.
#' @param intercept Optional `p x 1` matrix of initial intercept estimates.
#'   If `NULL` (default), computed alongside `betaEst`.
#' @param scaleVar Optional numeric vector of length `p` containing initial
#'   scale estimates. If `NULL` (default), computed internally.
#' @param xhat Optional `n x p` matrix of expected (cleaned) predictor values.
#'   If `NULL` (default), computed via a robust univariate regression step.
#' @param xtilde Optional `n x p` matrix of pre-processed predictor values.
#'   If `NULL` (default), derived from `xhat`.
#' @param maxituniv Maximum number of IRLS iterations per univariate
#'   sub-regression. Default is `1` (one-step re-weighting).
#' @param maxitscale Maximum number of fixed-point iterations for the M-scale
#'   estimate. Default is `100`.
#' @param wvalue Cut-off value used for cellwise outlier flagging. Default
#'   is `3`.
#' @param shoot_order Order in which predictors are cycled. `"default"` cycles
#'   in column order `1, ..., p`; `"robcor"` orders by decreasing robust
#'   (Kendall) correlation with `y`.
#' @param nlambda Number of values in the automatically constructed lambda
#'   grid. Default is `100`. Ignored when `lambda_grid` is supplied.
#' @param post Logical. If `TRUE` (default), the post-lasso approach is used
#'   (coefficients are not shrunk after support identification).
#' @param standardize Logical. If `TRUE` (default), predictors are robustly
#'   standardised before fitting and coefficients are back-transformed on
#'   output.
#' @param lambda_grid Optional numeric vector of sparsity parameter values to
#'   evaluate. If `NULL` (default), a grid is constructed automatically.
#' @param predset Optional integer vector specifying the initial predictor
#'   subset used for the MM starting-value fit.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{`coef`}{Numeric vector of length `p + 1`: intercept followed by
#'       `p` slope coefficients, selected by BIC using the M-scale.}
#'     \item{`coef_ln`}{Numeric vector of length `p + 1`: alternative
#'       coefficient estimates selected by BIC using the log M-scale.}
#'     \item{`weights`}{`n x p` matrix of robustness weights at the optimal
#'       lambda.}
#'     \item{`weights_ln`}{`n x p` robustness weights at the log-BIC optimal
#'       lambda.}
#'     \item{`lambda_opt`}{Optimal sparsity parameter selected by BIC.}
#'     \item{`lambdas`}{Full lambda grid that was evaluated.}
#'     \item{`iter`}{Number of outer shooting loop iterations at the optimal
#'       lambda.}
#'     \item{`betahats`}{List of slope coefficient vectors, one per lambda.}
#'     \item{`alphahats`}{List of intercept values, one per lambda.}
#'     \item{`fits`}{Full output of the internal model selection step,
#'       containing results for all lambda values.}
#'     \item{`start_flag`}{Logical matrix indicating which predictor cells
#'       were flagged as outliers during initialisation.}
#'     \item{`startfit`}{Output of the MM initialisation step.}
#'     \item{`xhat`}{Cleaned predictor matrix used internally.}
#'     \item{`xtilde`}{Pre-processed predictor matrix used internally.}
#'   }
#'
#' @references Bottmer, L., Croux, C., and Wilms, I. (2022). Sparse regression
#'   for large data sets with outliers. *Econometrics and Statistics*.
#'
#' @seealso [shooting()] for the non-sparse variant.
#'
#' @examples
#' set.seed(1)
#' n <- 50; p <- 20
#' X <- matrix(rnorm(n * p), n, p)
#' y <- X[, 1] + X[, 2] + rnorm(n)
#' fit <- sparseshooting(x = X, y = y, nlambda = 20)
#' fit$coef
#'
#' @export
sparseshooting <- function(x, y, k = 3.420, maxIteration = 100, tol = 10^-2,
  betaEst = NULL, intercept = NULL, scaleVar = NULL, xhat = NULL,
  xtilde = NULL, maxituniv = 1, maxitscale = 100, wvalue = 3,
  shoot_order = "default", nlambda = 100, post = TRUE, standardize = TRUE,
  lambda_grid = NULL, predset = NULL) {
  n <- nrow(x)
  p <- ncol(x)
  delta <- (1 - 3 / k^2 + 5 / k^4 - k^2 / 3) * pnorm(k) +
    (4 / (3 * k) - k / 3 - 5 / k^3) * dnorm(k) -
    1 / 2 + 3 / (2 * k^2) - 5 / (2 * k^4) + k^2 / 3

  if (standardize) {
    Xest <- Xestimfast(x, value = wvalue)
    Xclean <- apply(x, 2, Xinitftc, Xmatrix = x, Xest = Xest, value = wvalue)
    sx <- apply(Xclean, 2, sd)
    mx <- apply(Xclean, 2, mean)
    x <- scale(x, center = mx, scale = sx)
  }

  if (is.null(xhat)) {
    xhat <- Xestimfast(Xmatrix = x, value = wvalue)
  }
  if (is.null(xtilde)) {
    xtilde <- apply(X = x, 2, Xinitftc, Xmatrix = x, Xest = xhat)
  }
  xhat0 <- xhat
  xtilde0 <- xtilde
  start_flag <- xtilde != x

  if (shoot_order == "default") {
    order_variables <- 1:p
  }

  robcor_fit <- NULL
  if (shoot_order == "robcor") {
    robcor_fit <- robcorr(X = x, Y = y)
    order_variables <- robcor_fit$predcor
  }

  startfit <- NULL
  if (is.null(betaEst)) {
    startfit <- startvalue_MM(
      X = x, Y = y, value = wvalue, robcor_fit = robcor_fit,
      Xinit = xtilde, predset = predset
    )
    betaEst <- as.matrix(startfit$betaEst)
    intercept <- as.matrix(startfit$intercept)
    scaleVar <- startfit$scaleVar
  }

  betaEst <- as.matrix(betaEst[order_variables, 1])
  intercept <- as.matrix(intercept[order_variables, 1])
  scaleVar <- scaleVar[order_variables]
  xtilde <- xtilde[, order_variables]
  ytilde <- y - xtilde %*% betaEst
  x <- x[, order_variables]
  xhat <- xhat[, order_variables]

  wt <- matrix(NA, ncol = p, nrow = n)
  y <- as.matrix(y)

  if (is.null(lambda_grid)) {
    lambda_max <- get_lambda_max(
      xtilde, y, ytilde, x, xhat, betaEst, intercept, scaleVar,
      k, delta, maxIteration, tol, maxitscale, maxituniv, post,
      value = wvalue
    )
    lambda_grid <- getLambdaGrid(lambda_max, n, p, nlambda)
  }

  fit_regression <- selectModel(
    lambda_grid = lambda_grid, ytilde = ytilde, xtilde = xtilde, xhat = xhat,
    x = x, wt = wt, betaEst = betaEst, intercept = intercept,
    scaleVar = scaleVar, k = k, delta = delta, maxIteration = maxIteration,
    y = y, maxitscale = maxitscale, maxituniv = maxituniv, tol = tol,
    post = post, value = wvalue
  )

  betahat <- betahat_ln <- rep(NA, p)
  betahat[order_variables] <- fit_regression$coef_opt
  betahat_ln[order_variables] <- fit_regression$coef_opt_ln
  betahats <- lapply(
    fit_regression$get_coefs$betaEst,
    p = p, order_variables = order_variables,
    function(X, p, order_variables) {
      get_beta <- rep(NA, p)
      get_beta[order_variables] <- X
      return(get_beta)
    }
  )

  alphahat <- fit_regression$alpha_opt
  alphahat_ln <- fit_regression$alpha_opt_ln
  alphahats <- fit_regression$get_coefs$alpha
  weights <- weights_ln <- matrix(NA, n, p)
  weights[, order_variables] <- fit_regression$weights_opt
  weights_ln[, order_variables] <- fit_regression$weights_opt_ln

  if (standardize) {
    betahat <- betahat / sx
    betahat_ln <- betahat_ln / sx
    alphahat <- alphahat - sum(betahat * mx)
    alphahat_ln <- alphahat_ln - sum(betahat_ln * mx)

    betahats <- lapply(betahats, sx = sx, function(X, sx) X / sx)
    alphahats <- lapply(
      fit_regression$lambdas,
      lambdas = fit_regression$lambdas,
      alphas = alphahats,
      betas = betahats,
      mx = mx,
      function(X, lambdas, alphas, betas, mx) {
        index <- which(lambdas == X)
        alphas[[index]] - sum(betas[[index]] * mx)
      }
    )
  }

  list(
    "coef" = c(alphahat, betahat),
    "coef_ln" = c(alphahat_ln, betahat_ln),
    "weights" = weights,
    "weights_ln" = weights_ln,
    "lambda_opt" = fit_regression$lambda_opt,
    "lambdas" = fit_regression$lambdas,
    "iter" = fit_regression$iter_opt,
    "betahats" = betahats,
    "alphahats" = alphahats,
    "fits" = fit_regression,
    "start_flag" = start_flag,
    "startfit" = startfit,
    "xhat" = xhat0,
    "xtilde" = xtilde0
  )
}
