set.seed(9341)
n <- 40
p <- 10
X <- matrix(rnorm(n * p), n, p)
y <- X[, 1] + X[, 2] + rnorm(n)

test_that("sparseshooting returns expected structure", {
  fit <- sparseshooting(x = X, y = y, nlambda = 10)

  expect_type(fit, "list")
  expect_named(fit, c(
    "coef", "coef_ln", "weights", "weights_ln", "lambda_opt", "lambdas",
    "iter", "betahats", "alphahats", "fits", "start_flag", "startfit",
    "xhat", "xtilde"
  ))
  # intercept + p slopes
  expect_length(fit$coef, p + 1)
  expect_length(fit$coef_ln, p + 1)
  expect_equal(dim(fit$weights), c(n, p))
  expect_length(fit$lambdas, 11) # nlambda + 1 (p < n adds lambda=0)
  expect_true(is.numeric(fit$lambda_opt))
  expect_true(fit$lambda_opt >= 0)
})

test_that("sparseshooting coefficients are finite", {
  fit <- sparseshooting(x = X, y = y, nlambda = 10)
  expect_true(all(is.finite(fit$coef)))
})

test_that("sparseshooting recovers sparse signal", {
  fit <- sparseshooting(x = X, y = y, nlambda = 20)
  # First two predictors should have non-zero coefficients
  expect_gt(abs(fit$coef[2]), 0.1)
  expect_gt(abs(fit$coef[3]), 0.1)
})

test_that("shooting returns expected structure", {
  fit <- shooting(x = X, y = y)

  expect_type(fit, "list")
  expect_named(fit, c(
    "coef", "weights", "iter", "start_flag", "startfit", "xhat", "xtilde"
  ))
  expect_length(fit$coef, p + 1)
  expect_equal(dim(fit$weights), c(n, p))
  expect_true(all(is.finite(fit$coef)))
})

test_that("shooting recovers signal", {
  fit <- shooting(x = X, y = y)
  expect_gt(abs(fit$coef[2]), 0.3)
  expect_gt(abs(fit$coef[3]), 0.3)
})

test_that("standardize = FALSE runs without error", {
  expect_no_error(shooting(x = X, y = y, standardize = FALSE))
  expect_no_error(sparseshooting(x = X, y = y, nlambda = 5, standardize = FALSE))
})

test_that("p = n/2 no longer triggers a warning after kpred fix", {
  # Previously kpred = round(n/2) caused lmrob to fail at the breakdown boundary.
  # Fixed to kpred = n %/% 2 - 1, so p = n/2 now works cleanly.
  set.seed(1)
  X_big <- matrix(rnorm(100 * 50), 100, 50)
  y_big <- X_big[, 1] + rnorm(100)
  expect_no_warning(sparseshooting(x = X_big, y = y_big, nlambda = 5))
})
