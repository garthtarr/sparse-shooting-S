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

# Simulate data
set.seed(1)
n <- 100; p <- 50
X <- matrix(rnorm(n * p), n, p)
beta <- c(rep(1, 5), rep(0, 45))
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
- `robStandardize()` was incorrectly assumed to be from `robustbase`; it is in fact from `robustHD`. The AI initially replaced it with `median()` / `mad()`, but this was corrected after user feedback. `robustHD` has been added as an explicit dependency and `robStandardize()` reinstated.
- `stats` functions (`dnorm`, `mad`, `median`, `pnorm`, `sd`) added as explicit imports.
- `testthat` added to `Suggests`; `R (>= 4.1.0)` added to `Depends` (required for `\()` lambda syntax used in utils).

---

## Original authors

Lea Bottmer, Christophe Croux, Ines Wilms
