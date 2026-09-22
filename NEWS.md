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
