# Internal auxiliary functions for sparseshootS
# None of these are exported; they support sparseshooting() and shooting().

robcorr <- function(X, Y) {
  robcor <- apply(X, 2, pcaPP::cor.fk, y = Y)
  predcor <- order(abs(robcor), decreasing = TRUE)
  list("robcor" = robcor, "predcor" = predcor)
}

Xestimfast <- function(Xmatrix, value = 3) {
  cork <- abs(pcaPP::cor.fk(Xmatrix))
  diag(cork) <- 0
  cormax <- apply(cork, 2, which.max)

  MMuniv <- function(U, indexmax, Xdata) {
    i.variable <- which(apply(U == Xdata, 2, all) == TRUE)
    i.pred <- indexmax[i.variable]
    Xstand <- robustHD::robStandardize(U)
    mx <- attr(Xstand, "center")
    sx <- attr(Xstand, "scale")
    fit <- suppressWarnings(
      robustbase::lmrob(U ~ Xdata[, i.pred], setting = "KS2011")
    )
    Xstar <- fit$fitted.values
    xstars <- (Xstar - mx) / sx
    flag <- abs(xstars) > value
    Xstar * (1 - flag) + median(U) * flag
  }

  apply(Xmatrix, 2, MMuniv, indexmax = cormax, Xdata = Xmatrix)
}

Xinitftc <- function(U, Xmatrix, Xest, value = 3) {
  i.variable <- which(apply(U == Xmatrix, 2, all) == TRUE)
  xx <- robustHD::robStandardize(U)
  mx <- attr(xx, "center")
  sx <- attr(xx, "scale")
  critvalue <- value * sx
  flagX <- 1 * (abs(U - mx) > critvalue)
  (1 - flagX) * U + flagX * (Xest[, i.variable])
}

get_lambda_max <- function(xtilde, y, ytilde, x, xhat, betaEst, intercept,
  scaleVar, k, delta, maxIteration, tol, maxitscale, maxituniv, post,
  value) {
  cork <- apply(xtilde, 2, pcaPP::cor.fk, y = y)
  max_idx <- which.max(abs(cork))
  lambda_init <- 30 * abs(cork[max_idx])
  n <- nrow(x)
  p <- ncol(x)

  lambda_max <- lambda_init
  lambda_last_max <- lambda_init
  last_frac <- 0.5
  count <- 1

  for (i in 1:10) {
    wt <- matrix(NA, ncol = p, nrow = n)
    ft <- shootloop_sparse(
      ytilde = ytilde, xtilde = xtilde, x = x, wt = wt, Xexp = xhat,
      betaEst = betaEst, intercept = intercept, scaleVar = scaleVar,
      scaleVarOLD = scaleVar, k = k, delta = delta,
      maxIteration = maxIteration, tol = mad(y) * 10^-6,
      tolout = tol, tolscale = 10^-5, maxitscale = maxitscale, y = y,
      maxituniv = maxituniv, lambda = lambda_max, post = post, wcut = value
    )
    if (sum(ft$betaEst == 0) == p) {
      lambda_last_max <- lambda_max
      count <- 1
      last_frac <- 0.5
      lambda_max <- 0.5 * lambda_max
    } else {
      if (i == 1) {
        lambda_max <- lambda_last_max <- 10 * lambda_max
      } else {
        count <- count + 1
        last_frac <- last_frac + 1 / 2^count
        lambda_max <- last_frac * lambda_last_max
        if (last_frac >= 1) break
      }
    }
  }
  lambda_last_max
}

getLambdaGrid <- function(lambda_max, n, p, size) {
  lmin <- lambda_max / 10^5
  grid <- exp(seq(log(lambda_max), log(lmin), length = size))
  if (p < n) grid <- c(grid, 0)
  grid
}

selectModel <- function(lambda_grid, ytilde, xtilde, xhat, x, wt, betaEst,
  intercept, scaleVar, k, delta, maxIteration, y, maxitscale, maxituniv,
  tol, post, value) {
  n <- length(ytilde)
  p <- ncol(xtilde)

  get_coefs <- shootloop_sparse_lambdas(
    ytilde = ytilde, xtilde = xtilde, x = x, wt = wt, Xexp = xhat,
    betaEst = betaEst, intercept = intercept, scaleVar = scaleVar,
    scaleVarOLD = scaleVar, k = k, delta = delta,
    maxIteration = maxIteration, tol = mad(y) * 10^-6,
    tolout = tol, tolscale = 10^-5, maxitscale = maxitscale, y = y,
    maxituniv = maxituniv, lambda = lambda_grid, post = post, wcut = value
  )

  dfs <- unlist(lapply(get_coefs$betaEst, \(U) length(which(U != 0))))
  sigmahat <- unlist(get_coefs$sigmahat)
  bic_sigmas <- sigmahat^2 + dfs * (log(n) / n)
  bic_ln_sigmas <- log(sigmahat^2) + dfs * (log(n) / n)
  minimum_sigmas <- which.min(bic_sigmas)
  minimum_ln_sigmas <- which.min(bic_ln_sigmas)

  list(
    "get_coefs" = get_coefs,
    "lambdas" = lambda_grid,
    "coef_opt" = get_coefs$betaEst[[minimum_sigmas]],
    "lambda_opt" = lambda_grid[minimum_sigmas],
    "BIC_opt" = bic_sigmas[minimum_sigmas],
    "weights_opt" = get_coefs$weights[[minimum_sigmas]],
    "alpha_opt" = get_coefs$alpha[[minimum_sigmas]],
    "iter_opt" = get_coefs$iter[[minimum_sigmas]],
    "coef_opt_ln" = get_coefs$betaEst[[minimum_ln_sigmas]],
    "lambda_opt_ln" = lambda_grid[minimum_ln_sigmas],
    "BIC_opt_ln" = bic_sigmas[minimum_ln_sigmas],
    "weights_opt_ln" = get_coefs$weights[[minimum_ln_sigmas]],
    "alpha_opt_ln" = get_coefs$alpha[[minimum_ln_sigmas]],
    "iter_opt_ln" = get_coefs$iter[[minimum_ln_sigmas]]
  )
}

startvalue_MM <- function(X, Y, value, robcor_fit = NULL, Xinit = NULL,
  predset = NULL) {
  p <- ncol(X)
  n <- nrow(X)
  betaEst <- rep(0, p)

  if (is.null(predset)) {
    kpred <- min(n %/% 2 - 1, p)
    if (is.null(robcor_fit)) {
      robcor_fit <- robcorr(X = X, Y = Y)
    }
    predset <- robcor_fit$predcor[1:kpred]
  }

  fitinit <- suppressWarnings(
    robustbase::lmrob(
      Y ~ ., data = cbind(Y, data.frame(Xinit[, predset])),
      setting = "KS2011"
    )
  )

  if (!isTRUE(fitinit$converged) || fitinit$scale == 0) {
    warning(
      "The MM initialisation (lmrob) did not converge and returned a scale of ",
      "zero. This will cause all Tukey biweight weights in the shooting loop to ",
      "collapse to zero, so the penalisation is never applied and coefficients ",
      "will not be sparse.",
      call. = FALSE
    )
  }

  betaEst[predset] <- fitinit$coef[-1]
  intercept <- rep.int(fitinit$coef[1], p)
  scaleVar <- rep.int(fitinit$scale, p)

  list(
    "betaEst" = as.matrix(betaEst),
    "intercept" = as.matrix(intercept),
    "scaleVar" = scaleVar,
    "Xinit" = Xinit
  )
}
