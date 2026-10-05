# sempredict (development version)

* `cv_predict()` no longer stops with "subscript out of bounds" when there is a
  single outcome and the `"lm"` benchmark (the default) is used: with one
  response `lm()` is not a multivariate fit and `predict()` returns an unnamed
  vector, which is now reshaped to a named matrix.
* `plot_dispersion()` no longer errors when a column has no predictions in the
  drawn repetition (for example a PLSc `implied` column that was inadmissible
  on every fold); the column is dropped with a warning.
* `assess()` identifies the mean benchmark by its role, not its name, so a
  model column labelled `"mean"` (with the mean benchmark switched off) now
  gets a dispersion slope.
* `predict_oos()` and `cv_predict()` now refuse indicators listed in both
  `ynames` and `xnames`; an outcome used as its own predictor gave perfect,
  meaningless predictions.
* Tests: the lavaan oracle test no longer fails on lavaan 0.6-17, whose
  `as.matrix()` keeps the `"lavaan.matrix"` class.

# sempredict 0.1.3

* `sem_params()` for seminr PLSc models now aligns the implied construct
  correlations with the loadings by name. Version 0.1.2 multiplied them by
  position, so the model-implied rule (`predict_oos(construction = "implied")`)
  was wrong whenever the order in which `relationships()` mentions the
  constructs differed from their causal order. Models declared in causal order
  (all examples, the vignette and the tutorial) were not affected.
* `sem_params()` no longer errors ("missing value where TRUE/FALSE needed") when a
  PLSc solution has rho_A <= 0: the resulting non-finite loadings and paths are now
  reported as an inadmissible solution (`admissible = FALSE`, `reason` naming
  rho_A), as documented. Found in a simulation at n = 100 with weak structural paths.
* `predict_oos(construction = "chain")` now stops with a clear message instead of
  returning NaN predictions when the fitted weights, paths or loadings are not
  finite, so `cv_predict()` counts the fold as a failure.

# sempredict 0.1.2

* `cvpat()` now defaults to the one-sided test of Liengaard et al. (2021) and
  Sharma et al. (2023), `alternative = "greater"`: the alternative hypothesis is
  that the model's loss is lower than the benchmark's. Version 0.1.1 defaulted
  to a two-sided test, which halves nothing when the loss difference has the
  hypothesised sign but reports `p/2` instead of `1 - p/2` when it does not.
  Pass `alternative = "two.sided"` for the earlier behaviour.

# sempredict 0.1.1

* `sem_params()` for consistent PLS: a standardised loading above one (a
  negative residual variance) now makes the solution inadmissible, and `reason`
  names the items. Version 0.1.0 floored the residual variance at 1e-6, which
  turns the affected item into an error-free measure of its construct and lets
  it dominate the model-implied regression. On the tutorial's PoliticalDemocracy
  folds the model-implied regression is now unavailable on 82 of 200 training
  fits (0.1.0: 17); the ECSI example is unaffected. Other inadmissibility
  reasons note when a loading also exceeds one, so `failures()` reproduces the
  tutorial's failure table.

# sempredict 0.1.0

* First prototype, extracted from the companion code of the Danks & Sharma
  Teacher's Corner tutorial. `sem_params()`, `predict_oos()`, `sem_method()`,
  `cv_predict()`, `assess()`, `failures()`, `cvpat()`, `cvpat_table()`,
  `plot_dispersion()`.
