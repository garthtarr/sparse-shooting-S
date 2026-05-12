#' Shooting S-Estimator
#'
#' Computes the (non-sparse) shooting S-estimator for robust regression in
#' high-dimensional settings with cellwise contamination.
#'
#' @inheritParams sparseshooting
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{`coef`}{Numeric vector of length `p + 1`: intercept followed by
#'       `p` slope coefficients.}
#'     \item{`weights`}{`n x p` matrix of robustness weights.}
#'     \item{`iter`}{Number of outer shooting loop iterations.}
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
#' @seealso [sparseshooting()] for the sparse variant with automatic lambda
#'   selection.
#'
#' @examples
#' set.seed(1)
#' n <- 50; p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' y <- X[, 1] + X[, 2] + rnorm(n)
#' fit <- shooting(x = X, y = y)
#' fit$coef
#'
#' @export
shooting <- function(x, y, k = 3.420, maxIteration = 100, tol = 10^-2,
  betaEst = NULL, intercept = NULL, scaleVar = NULL, xhat = NULL,
  xtilde = NULL, maxituniv = 1, maxitscale = 100, wvalue = 3,
  shoot_order = "default", standardize = TRUE, predset = NULL) {
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

  shootfit <- shootloop(
    ytilde = ytilde, xtilde = xtilde, x = x, wt = wt, Xexp = xhat,
    betaEst = betaEst, intercept = intercept, scaleVar = scaleVar,
    scaleVarOLD = scaleVar, k = k, delta = delta,
    maxIteration = maxIteration, tol = mad(y) * 10^-6,
    tolout = tol, tolscale = 10^-5, maxitscale = maxitscale,
    y = y, maxituniv = 1, wcut = wvalue
  )

  betahat <- rep(NA, p)
  betahat[order_variables] <- shootfit$betaEst
  alphahat <- shootfit$alpha
  weights <- matrix(NA, n, p)
  weights[, order_variables] <- shootfit$weights

  if (standardize) {
    betahat <- betahat / sx
    alphahat <- alphahat - sum(betahat * mx)
  }

  list(
    "coef" = c(alphahat, betahat),
    "weights" = weights,
    "iter" = shootfit$iter,
    "start_flag" = start_flag,
    "startfit" = startfit,
    "xhat" = xhat0,
    "xtilde" = xtilde0
  )
}
