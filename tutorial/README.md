# Companion code — "Assessing Out-of-Sample Prediction in CB-SEM, PLS-SEM, and PLSc with lavaan and seminr"
[Authors to be added; author order not yet decided]. Teacher's Corner, *Structural Equation Modeling* (submitted 2026).

Everything the article reports regenerates from this folder with seed 20260825.
Run from this folder's root (the scripts write a relative `results/tutorial/`):

```
Rscript R/tutorial/00_quickstart.R            # compact ML / PLS / PLSc example (Step 3-4 chunks)
Rscript R/tutorial/01_political_democracy.R   # Table 2, the failure counts and the residual-covariance check (~40 s)
Rscript R/tutorial/02_mobi_low_reliability.R  # Table 3 (~30 s)
Rscript R/tutorial/03_slope_finite_sample.R   # Supplement: finite-sample dispersion-slope reference (~2 min)
Rscript R/tutorial/04_cvpat.R                 # CVPAT on one 10-fold repetition (reads the .rds from 01 and 02)
Rscript R/tutorial/05_figures.R               # article Figures 1 and 2 (PNG + 300-dpi TIFF in figures/; reads the .rds from 02)
```

`R/tutorial/00_helpers.R` holds the shared machinery: the construct-score chain, the
covariance-based best linear predictor, PLSc parameter extraction with the admissibility check,
repeated k-fold cross-validation with identical folds for every method, RMSE with Monte Carlo
standard errors, overfit ratio, dispersion slope, and common-fold rescoring.

Data: `lavaan::PoliticalDemocracy` and `seminr::mobi` (shipped with the packages).
Software: R >= 4.5, lavaan >= 0.6-13 (for `lavPredictY()`), seminr 2.6.0 (the CRAN release, 2.5.0, reproduces all outputs; verified 22 Sep 2026), MASS.
`results/sessionInfo.txt` records the versions used for the reported numbers; the text outputs in
`results/tutorial/` are the runs quoted in the article.

The supporting 128-cell simulation (supplement section S4) is archived separately in the authors'
replication package for the simulation study; the simulation numbers quoted in the article regenerate from `R/analysis/s4_tutorial_simulation_summary.R` in that package, not from this folder.
