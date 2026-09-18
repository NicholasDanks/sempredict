test_that("lavaan params carry the implied moments and are admissible", {
  fit <- lavaan::sem(pd_syntax, data = pd, meanstructure = TRUE)
  p <- sem_params(fit)
  expect_s3_class(p, "sem_params")
  expect_true(p$admissible)
  expect_null(p$chain)
  expect_equal(dim(p$Sigma), c(11, 11))
  expect_equal(names(p$mu), colnames(p$Sigma))
  imp <- lavaan::lavInspect(fit, "implied")
  expect_equal(unname(p$mu), as.numeric(imp$mean))
})

test_that("lavaan params refuse a model without a mean structure", {
  fit <- lavaan::sem(pd_syntax, data = pd)
  expect_error(sem_params(fit), "meanstructure")
})

test_that("seminr PLS params: composite metric, rho_A = 1, chain available", {
  m <- quiet_pls(pd, pd_mm_pls, pd_sm)
  p <- sem_params(m)
  expect_equal(p$estimator, "PLS")
  expect_true(p$admissible)
  expect_equal(unname(p$rhoA), c(1, 1, 1))
  expect_equal(p$chain$order, c("ind60", "dem60", "dem65"))
  expect_equal(names(p$chain$blocks$dem65), NULL)
  expect_equal(p$chain$blocks$dem65, paste0("y", 5:8))
  # implied Sigma is on the raw scale: its diagonal is close to the sample variances
  # scaled by the model (unit variance in z metric times sd^2 when Theta fills the gap)
  expect_equal(unname(diag(p$Sigma)), unname(p$sd^2), tolerance = 1e-6)
})

test_that("seminr PLSc params: rho_A in (0,1], corrected correlations, admissibility", {
  m <- quiet_pls(pd, pd_mm_plsc, pd_sm)
  p <- sem_params(m)
  expect_equal(p$estimator, "PLSc")
  expect_true(p$admissible)
  expect_true(all(p$rhoA > 0 & p$rhoA <= 1))
  # the tutorial reports rho_A(ind60) = .954 on the full sample
  expect_equal(unname(p$rhoA["ind60"]), 0.954, tolerance = 1e-3)
  # Phi is a correlation matrix consistent with the structural model
  expect_equal(unname(diag(p$Phi)), c(1, 1, 1))
  expect_true(min(eigen(p$Phi, only.values = TRUE)$values) > 0)
})

test_that("an inadmissible disattenuated correlation is reported, not raised", {
  # Build a two-construct data set whose composites correlate .95 while the
  # blocks are unreliable, so the disattenuated correlation exceeds one.
  set.seed(1)
  n <- 200
  f <- rnorm(n)
  d <- data.frame(a1 = f + rnorm(n, sd = 1.2), a2 = f + rnorm(n, sd = 1.2), a3 = f + rnorm(n, sd = 1.2))
  s <- rowMeans(d)
  d$b1 <- s + rnorm(n, sd = .3); d$b2 <- s + rnorm(n, sd = .3); d$b3 <- s + rnorm(n, sd = .3)
  mm <- seminr::constructs(seminr::reflective("A", paste0("a", 1:3)), seminr::reflective("B", paste0("b", 1:3)))
  sm <- seminr::relationships(seminr::paths(from = "A", to = "B"))
  m <- quiet_pls(d, mm, sm)
  p <- sem_params(m)
  expect_false(p$admissible)
  expect_match(p$reason, "inadmissible")
  expect_null(p$Sigma)
  expect_false(is.null(p$chain))
})

test_that("topological order handles a non-model order and rejects cycles", {
  B <- matrix(0, 3, 3, dimnames = list(c("c", "a", "b"), c("c", "a", "b")))
  B["a", "b"] <- .5; B["b", "c"] <- .4
  expect_equal(sempredict:::topo_order(B), c("a", "b", "c"))
  B["c", "a"] <- .1
  expect_error(sempredict:::topo_order(B), "recursive")
})

# The tutorial's fold stream: seed 20260825, one draw of fold labels per repetition.
tutorial_folds <- function(n, reps = 20, k = 10, seed = 20260825) {
  set.seed(seed)
  lapply(seq_len(reps), function(r) sample(rep(seq_len(k), length.out = n)))
}

test_that("a PLSc loading above one is inadmissible, not floored", {
  # Repetition 1, first training fold whose PLSc solution has a standardised
  # loading above one (a negative residual variance).
  folds <- tutorial_folds(nrow(pd))[[1]]
  hit <- NULL
  for (k in 1:10) {
    m <- quiet_pls(pd[folds != k, ], pd_mm_plsc, pd_sm)
    if (any(rowSums(m$outer_loadings^2) > 1)) { hit <- m; break }
  }
  expect_false(is.null(hit))
  p <- sem_params(hit)
  expect_false(p$admissible)
  expect_match(p$reason, "loading above one")
  expect_match(p$reason, "x1|x2")           # the offending item is named
  expect_null(p$Sigma)
  expect_false(is.null(p$chain))            # the score chain stays available
  expect_error(predict_oos(p, newdata = pd[folds == k, ], xnames = pd_x_ea, ynames = pd_y,
                           construction = "implied"))
})

test_that("plain PLS is never flagged: composite loadings cannot exceed one", {
  p <- sem_params(quiet_pls(pd, pd_mm_pls, pd_sm))
  expect_true(p$admissible)
  expect_true(all(diag(p$Theta) >= 0))
})

test_that("tutorial oracle: 82 of 200 PoliticalDemocracy folds are inadmissible", {
  skip_on_cran()
  reasons <- character(0)
  for (folds in tutorial_folds(nrow(pd)))
    for (k in 1:10) {
      p <- sem_params(quiet_pls(pd[folds != k, ], pd_mm_plsc, pd_sm))
      if (!p$admissible) reasons <- c(reasons, p$reason)
    }
  expect_length(reasons, 82)
  expect_equal(sum(grepl("^inadmissible: standardised loading above one", reasons)), 65)
  expect_equal(sum(grepl("loading", reasons)), 73)   # 65 + 8 that also fail another check
  expect_equal(sum(grepl("not positive definite", reasons)), 15)
  expect_equal(sum(grepl("correlation >= 1", reasons)), 2)
})
