---
title: 'sempredict: Out-of-sample prediction assessment for structural equation models across estimators'
tags:
  - R
  - structural equation modeling
  - partial least squares
  - out-of-sample prediction
  - cross-validation
authors:
  # Author order not yet decided; listed alphabetically for now.
  - name: Nicholas P. Danks
    orcid: 0000-0001-6902-2708
    corresponding: true
    affiliation: 1
  - name: Pratyush N. Sharma
    orcid: 0000-0003-3913-0390
    affiliation: 2
affiliations:
  - name: Trinity Business School, Trinity College Dublin, Ireland
    index: 1
  - name: Culverhouse College of Business, The University of Alabama, USA
    index: 2
date: 17 September 2026
bibliography: paper.bib
---

# Summary

Structural equation models (SEM) are usually judged by how well they fit the
sample they were estimated on. A different question is how well a fitted model
predicts the outcomes of cases it has never seen. `sempredict` is an R package
that answers that question in one consistent way for the three estimators
applied researchers most often choose between: maximum-likelihood
covariance-based SEM fitted with `lavaan` [@rosseel2012], partial least squares
(PLS) fitted with `seminr` [@ray2024], and consistent PLS (PLSc; @dijkstra2015)
fitted with `seminr`'s reflective constructs. The package estimates nothing
itself. It turns any of these fits into case-level predictions of observed
outcome indicators from observed predictor indicators, cross-validates every
method on identical folds so that all comparisons are paired, and reports the
comparison with the diagnostics developed in the PLS prediction literature:
root mean squared error with a Monte Carlo standard error, benchmark
comparisons, the cross-validated predictive ability test (CVPAT;
@liengaard2021), an overfit ratio, a dispersion slope, and a log of every fold
on which a method failed and why.

Two ways of turning a fitted model into predictions are implemented. The
*construct-score chain* is the logic of PLSpredict [@shmueli2016;
@shmueli2019]: predictor items are weighted into construct scores, carried
through the path coefficients, and mapped back to the outcome items through the
loadings. The *model-implied regression* is the best linear prediction of the
outcome items from the predictor items under the fitted model's mean vector and
covariance matrix [@derooij2023], which for `lavaan` coincides with
`lavPredictY()`. For PLSc the package builds the implied covariances from the
corrected loadings and the reliability-corrected construct correlations and
checks the result for admissibility before predicting.

# Statement of need

Prediction has become a standard part of PLS-SEM reporting, with guidelines,
model-selection criteria and a paired test [@shmueli2019; @sharma2021;
@liengaard2021; @danks2024], while covariance-based SEM acquired a principled
out-of-sample predictor only recently [@derooij2023]. The two traditions have
developed in separate software with separate conventions, so a researcher who
wants to compare, on the same held-out cases, what an ML fit and a PLS fit of
the same model predict has had to write the harness by hand. `sempredict`
supplies that harness.

A second need is specific to consistent PLS. PLSc corrects loadings and path
coefficients for measurement unreliability but leaves the construct scores in
the composite metric. Applying the corrected coefficients to those scores, as
the construct-score chain does, over-disperses the predictions by roughly the
reciprocal square root of the antecedent's reliability, and on low-reliability
data the PLSc chain predicts worse than an ordinary regression on the raw
indicators. @cheah2026 note that PLSpredict and CVPAT "are not readily
applicable in PLSc-SEM". `sempredict` resolves this by predicting from the
PLSc solution through its implied covariances, in which scores and parameters
share a metric, and by reporting the solution as unavailable, rather than
silently wrong, whenever the reliability correction produces an inadmissible
construct correlation matrix. In the package's tests, on the ECSI customer
satisfaction data shipped with `seminr`, the chain's dispersion slope of
observed on predicted is 0.63 while the model-implied regression from the same
PLSc fits returns 0.92.

Existing software covers parts of this. `lavaan` predicts from one fitted
covariance model; `seminr::predict_pls()` implements PLSpredict for PLS;
`cSEM` [@rademaker2020] offers PLSpredict-style prediction with several
composite estimators and benchmarks. None runs several estimators through one
fold scheme with one set of diagnostics, and none provides the PLSc
model-implied regression with an admissibility check. `sempredict` is the
companion software for a tutorial on out-of-sample prediction across SEM
estimators [@danks2026] and reproduces every number in it.

# Design

The package has three layers, each a small set of functions.

1. **Extraction.** `sem_params()` is an S3 generic. For a `lavaan` fit it stores
   the model-implied mean vector and covariance matrix. For a `seminr` fit it
   stores the outer weights, path coefficients and loadings, derives the causal
   order of the constructs from the path matrix, and, when any construct is
   reflective, corrects the construct correlations by each construct's $\rho_A$,
   propagates them through the recursive structural model, and records whether
   the result is admissible.
2. **Prediction.** `predict_oos()` produces a matrix of predictions for new
   cases by either construction. The chain is available for `seminr` fits; the
   model-implied regression is available for any fit whose parameters are
   admissible.
3. **Assessment.** `sem_method()` wraps a fitting function with the
   constructions to apply; `cv_predict()` draws the fold assignments once, refits
   every method on every training fold, predicts the held-out fold, and stores
   predictions, squared errors, in-sample error and failures. `assess()`
   summarises RMSE, its Monte Carlo standard error across repetitions, the
   dispersion slope, the overfit ratio and failed folds, optionally rescoring all
   methods on the folds where a chosen method was available. `cvpat()` applies
   the paired test on one repetition; `failures()` lists reasons; and
   `plot_dispersion()` draws observed against predicted with the fitted slope.

```r
library(sempredict)
d  <- seminr::mobi
mm <- seminr::constructs(seminr::reflective("EXP", paste0("CUEX", 1:3)),
                         seminr::reflective("SAT", paste0("CUSA", 1:3)))
sm <- seminr::relationships(seminr::paths(from = "EXP", to = "SAT"))
plsc <- sem_method(function(tr) seminr::estimate_pls(tr, mm, sm),
                   construction = c("chain", "implied"))
cv <- cv_predict(list(PLSc = plsc), d, ynames = paste0("CUSA", 1:3),
                 xnames = paste0("CUEX", 1:3), k = 10, reps = 20, seed = 20260825)
assess(cv)
cvpat_table(cv, benchmark = "lm")
```

Correctness is tested against oracles rather than against the package's own
output: predictions from `lavaan` parameters equal `lavPredictY()` to
$10^{-10}$; construct scores reconstructed for the chain equal `seminr`'s own;
the ratio of the PLSc chain's dispersion to the model-implied regression's on
the ECSI data equals $1/\sqrt{\rho_A}$ to three decimals; and a twenty-times
ten-fold run with a fixed seed reproduces the tutorial's published table to
four decimals.

# Acknowledgements

The prediction constructions follow @derooij2023 and @shmueli2016; the
assessment diagnostics follow @shmueli2019, @liengaard2021 and @danks2024.

# References
