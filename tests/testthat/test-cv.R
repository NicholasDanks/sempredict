pls_method <- function(mm, sm, construction = "chain")
  sem_method(function(tr) quiet_pls(tr, mm, sm), construction = construction)
lav_method <- function(syntax)
  sem_method(function(tr) lavaan::sem(syntax, data = tr, meanstructure = TRUE, warn = FALSE))

test_that("cv_predict is reproducible from the seed and uses one fold set for every column", {
  ms <- list(ML = lav_method(pd_syntax), PLS = pls_method(pd_mm_pls, pd_sm))
  a <- cv_predict(ms, pd, pd_y, pd_x_da, k = 5, reps = 2, seed = 7)
  b <- cv_predict(ms, pd, pd_y, pd_x_da, k = 5, reps = 2, seed = 7)
  expect_equal(a$folds, b$folds)
  expect_equal(a$se_obs, b$se_obs)
  expect_equal(dim(a$folds), c(75, 2))
  expect_equal(a$columns, c("ML", "PLS", "mean", "lm"))
  # every case predicted once per repetition by every column
  expect_true(all(is.finite(a$se_obs)))
  expect_equal(unname(a$fail), c(0L, 0L, 0L, 0L))
})

test_that("a method with two constructions fits once and yields two columns", {
  ms <- list(PLSc = pls_method(pd_mm_plsc, pd_sm, construction = c("chain", "implied")))
  cv <- cv_predict(ms, pd, pd_y, pd_x_da, k = 5, reps = 1, seed = 3, benchmarks = NULL)
  expect_equal(cv$columns, c("PLSc_chain", "PLSc_implied"))
  expect_true(all(is.finite(cv$se_obs[, , "PLSc_chain"])))
})

test_that("failures are counted with reasons and do not stop the run", {
  broken <- sem_method(function(tr) stop("boom"))
  ms <- list(ML = lav_method(pd_syntax), bad = broken)
  cv <- cv_predict(ms, pd, pd_y, pd_x_da, k = 5, reps = 1, seed = 1)
  expect_equal(unname(cv$fail["bad"]), 5L)
  expect_true(all(is.na(cv$se_obs[, , "bad"])))
  f <- failures(cv)
  expect_equal(f$reason[f$method == "bad"], "boom")
  expect_equal(f$count[f$method == "bad"], 5L)
})

test_that("assess() returns the expected columns and common-fold rescoring restricts cells", {
  ms <- list(ML = lav_method(pd_syntax), bad = sem_method(function(tr) stop("boom")))
  cv <- cv_predict(ms, pd, pd_y, pd_x_da, k = 5, reps = 2, seed = 1)
  a <- assess(cv)
  expect_s3_class(a, "sem_assessment")
  expect_equal(names(a), c("method", "RMSE", "RMSE_MCSE", "dispersion_slope", "slope_MCSE", "OFR", "failed_folds"))
  expect_true(is.na(a$RMSE[a$method == "bad"]))
  expect_true(a$RMSE[a$method == "ML"] < a$RMSE[a$method == "mean"])
  expect_true(is.na(a$dispersion_slope[a$method == "mean"]))
  expect_equal(attr(a, "n_cells"), 150L)
  b <- assess(cv, common = "bad")
  expect_equal(attr(b, "n_cells"), 0L)
})

test_that("cvpat() reproduces a hand-computed paired t-test", {
  ms <- list(ML = lav_method(pd_syntax))
  cv <- cv_predict(ms, pd, pd_y, pd_x_da, k = 5, reps = 1, seed = 11)
  out <- cvpat(cv, "ML", "lm")
  d <- cv$se_obs[, 1, "lm"] - cv$se_obs[, 1, "ML"]
  tt <- t.test(d)
  expect_equal(out$d, mean(d))
  expect_equal(out$t, unname(tt$statistic))
  expect_equal(out$p, tt$p.value)
  expect_equal(out$n, 75L)
  tab <- cvpat_table(cv, "mean")
  expect_equal(tab$model, c("ML", "lm"))
})

test_that("tutorial oracle: mobi EXP -> SAT, 20 x 10-fold, seed 20260825 (Table 3 of the article)", {
  skip_on_cran()
  ms <- list(PLSc = pls_method(mobi_mm_plsc, mobi_sm, construction = c("chain", "implied")),
             PLS = pls_method(mobi_mm_pls, mobi_sm, construction = "chain"))
  cv <- cv_predict(ms, mobi, mobi_y, mobi_x, k = 10, reps = 20, seed = 20260825)
  a <- assess(cv)
  get <- function(col, v) a[[v]][a$method == col]
  expect_equal(get("PLSc_chain", "RMSE"), 1.4926, tolerance = 5e-4)
  expect_equal(get("PLSc_chain", "dispersion_slope"), 0.632, tolerance = 2e-3)
  expect_equal(get("PLSc_implied", "RMSE"), 1.4538, tolerance = 5e-4)
  expect_equal(get("PLSc_implied", "dispersion_slope"), 0.917, tolerance = 2e-3)
  expect_equal(get("PLS", "RMSE"), 1.4497, tolerance = 5e-4)
  expect_equal(get("lm", "RMSE"), 1.4610, tolerance = 5e-4)
  expect_equal(get("mean", "RMSE"), 1.5852, tolerance = 5e-4)
  expect_equal(get("PLSc_chain", "RMSE_MCSE"), 0.0012, tolerance = 2e-4)
  expect_equal(unname(cv$fail), rep(0L, 5))
  # CVPAT on repetition 1: the PLSc chain is the only row nominally worse than lm
  tab <- cvpat_table(cv, "lm", sensitivity = FALSE)
  expect_true(tab$d[tab$model == "PLSc_chain"] < 0)
  expect_equal(round(tab$p[tab$model == "PLSc_chain"], 3), 0.067)
  expect_equal(round(tab$p[tab$model == "PLSc_implied"], 3), 0.003)
})

test_that("plot_dispersion returns a ggplot", {
  skip_if_not_installed("ggplot2")
  ms <- list(PLS = pls_method(mobi_mm_pls, mobi_sm))
  cv <- cv_predict(ms, mobi, mobi_y, mobi_x, k = 5, reps = 1, seed = 2, benchmarks = NULL)
  expect_s3_class(plot_dispersion(cv), "ggplot")
})
