# sempredict

<!-- badges: start -->
[![R-CMD-check](https://github.com/NicholasDanks/sempredict/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/NicholasDanks/sempredict/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Out-of-sample prediction assessment for structural equation models, across
estimators.

Fit the model with **lavaan** (ML), **seminr** (PLS) or seminr with
`reflective()` constructs (consistent PLS, PLSc). sempredict turns any of these
fits into case-level predictions of observed outcome indicators from observed
predictor indicators, by either of two constructions, and cross-validates every
method on identical folds so that all comparisons are paired.

* **Model-implied regression**: the best linear prediction from the fitted
  model's mean vector and covariance matrix (de Rooij et al., 2023; `lavPredictY()`
  in lavaan). For PLSc the implied covariances are built from the corrected
  loadings and construct correlations and checked for admissibility.
* **Construct-score chain**: the PLSpredict logic (items → scores → paths →
  items). Under PLSc this applies factor-metric coefficients to
  composite-metric scores and over-disperses the predictions.

Reported per method: RMSE with a Monte Carlo standard error across
repetitions, the dispersion slope of observed on predicted, the overfit ratio,
mean and linear-regression benchmarks, the cross-validated predictive ability
test (CVPAT; Liengaard et al., 2021), and every prediction failure with its
reason.

## Installation

```r
# install.packages("remotes")
remotes::install_github("NicholasDanks/sempredict")
```

## Example

```r
library(sempredict)
d <- seminr::mobi
mm <- seminr::constructs(seminr::reflective("EXP", paste0("CUEX", 1:3)),
                         seminr::reflective("SAT", paste0("CUSA", 1:3)))
sm <- seminr::relationships(seminr::paths(from = "EXP", to = "SAT"))

plsc <- sem_method(function(tr) seminr::estimate_pls(tr, mm, sm),
                   construction = c("chain", "implied"))
cv <- cv_predict(list(PLSc = plsc), d,
                 ynames = paste0("CUSA", 1:3), xnames = paste0("CUEX", 1:3),
                 k = 10, reps = 20, seed = 20260825)
assess(cv)
cvpat_table(cv, benchmark = "lm")
plot_dispersion(cv)
```

See `vignette("sempredict")` for the full workflow on PoliticalDemocracy with
ML, PLS and PLSc side by side.

## Relation to other packages

lavaan provides `lavPredictY()` for a single fitted model; seminr provides
`predict_pls()` for PLS; cSEM provides `predict()` with several composite
estimators. sempredict does not estimate anything. It supplies the
cross-estimator harness with identical folds, the PLSc model-implied
regression with its admissibility check, and the assessment layer.

## Status

Prototype extracted from the companion code of Danks & Sharma, "Assessing
out-of-sample prediction in CB-SEM, PLS-SEM, and PLSc with lavaan and seminr"
(Teacher's Corner, *Structural Equation Modeling*, in preparation). The API may
change before a CRAN release.

## References

de Rooij, M., et al. (2023). SEM-based out-of-sample predictions. *Structural Equation Modeling*, 30(1), 132–148.
Liengaard, B. D., et al. (2021). Prediction: Coveted, yet forsaken? *Decision Sciences*, 52(2), 362–392.
Shmueli, G., et al. (2019). Predictive model assessment in PLS-SEM: Guidelines for using PLSpredict. *European Journal of Marketing*, 53(11), 2322–2347.
