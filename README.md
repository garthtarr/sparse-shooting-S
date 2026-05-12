# sparseshootS

R package implementing the **sparse shooting S-estimator** and **shooting S-estimator** for robust regression in high-dimensional settings with cellwise outliers.

This package is based on the method described in:

> Bottmer, L., Croux, C., and Wilms, I. (2022). Sparse regression for large data sets with outliers. *Econometrics and Statistics*.

---

## ⚠️ AI Disclaimer

**This R package was created by an AI assistant (Posit Assistant / Claude) without human involvement in the packaging process.**

The original scientific code was written by Lea Bottmer, Christophe Croux, and Ines Wilms. The AI was given the original source files (`sparseshootS.R`, `sparseshootS.cpp`, and `testscript.R`) and independently:

- Diagnosed and fixed compilation errors in the C++ source
- Restructured the code into a proper R package
- Wrote roxygen2 documentation
- Wrote unit tests
- Resolved dependency issues

**Any errors, omissions, or inaccuracies introduced during the packaging process are solely the responsibility of the AI and are not the fault of the original authors.** In particular:

- The documentation was written by the AI based on code inspection and may not accurately reflect the authors' intentions in all cases.
- Numerical behaviour of the original algorithms has not been independently verified beyond basic sanity checks.
- The unit tests were written by the AI and may not cover all edge cases.

If you find an error that is present in the original paper or original source files, please contact the original authors. If you find a packaging error (documentation, structure, build configuration), it was introduced by the AI.

---

## Installation

```r
pak::pak("garthtarr/sparse-shooting-S")
```

## Usage

```r
library(sparseshootS)

# Simulate data — keep p well below n/2
set.seed(1)
n <- 100
p <- 15
X <- matrix(rnorm(n * p), n, p)
beta <- c(rep(1, 3), rep(0, p-3))
y <- X %*% beta + rnorm(n)

# Sparse shooting S-estimator (automatic lambda selection via BIC)
fit_sparse <- sparseshooting(x = X, y = y, nlambda = 100)
fit_sparse$coef        # intercept + p slope coefficients
fit_sparse$lambda_opt  # selected sparsity parameter

# Non-sparse shooting S-estimator
fit <- shooting(x = X, y = y)
fit$coef
```

## What the functions do

### `sparseshooting()`

Fits a robust sparse linear regression model using the shooting S-estimator. The algorithm:

1. Robustly cleans the predictor matrix to handle cellwise outliers.
2. Initialises coefficients using an MM-estimator on the most correlated predictors.
3. Cycles through predictors using a shooting (coordinate descent) loop with IRLS and an M-scale objective.
4. Evaluates a grid of sparsity parameters (lambda) and selects the optimal value using a BIC criterion based on the residual M-scale.

### `shooting()`

The non-sparse variant: same shooting loop without the lasso-type penalisation.

---

## Known limitations

### Degenerate initialisation when `p ≥ n/2`

The algorithm initialises coefficients using an MM-estimator (`robustbase::lmrob`) fitted on the `kpred = min(round(n/2), p)` most robustly correlated predictors. When `p` is large relative to `n` — specifically when `kpred` approaches `n/2` — this MM fit can overfit the initialisation data, driving residuals to near zero and therefore the estimated scale (`lmrob$scale`) to zero.

A scale of zero propagates into the Tukey biweight weights used in the C++ shooting loop: every standardised residual becomes `Inf`, every weight becomes zero, the loop exits immediately via a `sumwjt == 0` guard, and `betaEst` is returned unchanged from the (degenerate) MM starting values. The penalisation is never applied, so **`sparseshooting()` returns all non-zero coefficients identical across the entire lambda grid** — as if no regularisation had been performed.

This was confirmed during the AI packaging process: with `n = 100` and `p = 50`, `lmrob` consistently returns `scale = 0` and `converged = FALSE`. With `p = 10` and `n = 100`, the algorithm works correctly and produces sparse solutions.

**Practical guidance:** ensure `p` is comfortably below `n/2`. The issue is in the initialisation heuristic in `startvalue_MM()`, not in the core penalisation logic (which was verified to work correctly in isolation). A future fix could cap `kpred` more conservatively or fall back to a simpler starting value when `lmrob` fails to converge.

A warning is raised automatically when this condition is detected:

```
Warning: The MM initialisation (lmrob) did not converge and returned a scale of
zero. This will cause all Tukey biweight weights in the shooting loop to collapse
to zero, so the penalisation is never applied and coefficients will not be sparse.
This typically occurs when p >= n/2. Try reducing p or increasing n.
```

---

## AI packaging process

The following steps were taken by the AI to convert the original scripts into this package:

### Bugs fixed in `sparseshootS.cpp`

1. **Typo in Rcpp export attribute** (`// [[Rcpp::export}]]` → `// [[Rcpp::export]]`) on the `sign` function.
2. **Wrong C++ standard plugin** — `// [[Rcpp::plugins(cpp11)]]` was forcing C++11 mode, but the installed RcppArmadillo requires C++14+. Removing the plugin lets `sourceCpp` use the system default (C++17).
3. **Unqualified `accu` calls** — `accu(wjt)` was replaced with `arma::accu(wjt)` in two places, as the `using namespace arma` directive is silently ignored (noted in the original code's own comments) and newer Armadillo no longer resolved the bare name.
4. **Extraneous parentheses** in two `if` conditions (compiler warnings, not errors).

### Packaging changes

- R code split into `R/sparseshooting.R`, `R/shooting.R`, `R/utils.R`, and `R/sparseshootS-package.R`.
- C++ source moved to `src/`.
- `Rcpp::compileAttributes()` run to generate `R/RcppExports.R` and `src/RcppExports.cpp`.
- `stats` functions (`dnorm`, `mad`, `median`, `pnorm`, `sd`) added as explicit imports.
- `testthat` added to `Suggests`; `R (>= 4.1.0)` added to `Depends` (required for `\()` lambda syntax used in utils).

---

## Original authors

Lea Bottmer, Christophe Croux, Ines Wilms
