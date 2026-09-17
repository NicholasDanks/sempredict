test_that("model-implied regression from lavaan equals lavPredictY (oracle)", {
  set.seed(20260825)
  id <- sample(seq_len(nrow(pd)), 60)
  fit <- lavaan::sem(pd_syntax, data = pd[id, ], meanstructure = TRUE)
  p <- sem_params(fit)
  for (x in list(pd_x_da, pd_x_ea)) {
    ours <- predict_oos(p, pd[-id, ], pd_y, x, construction = "implied")
    theirs <- lavaan::lavPredictY(fit, newdata = pd[-id, ], ynames = pd_y, xnames = x)
    expect_equal(unname(ours), unname(as.matrix(theirs)), tolerance = 1e-10)
  }
})

test_that("the chain is unavailable for lavaan and errors informatively", {
  fit <- lavaan::sem(pd_syntax, data = pd, meanstructure = TRUE)
  expect_error(predict_oos(sem_params(fit), pd, pd_y, pd_x_da, construction = "chain"), "seminr")
})

test_that("seminr chain: standardised construct scores reproduce seminr's own", {
  m <- quiet_pls(pd, pd_mm_pls, pd_sm)
  p <- sem_params(m)
  # Scores of an exogenous construct scored from its own items on the training
  # data must equal seminr's construct scores (both standardised on the sample).
  Z <- scale(as.matrix(pd[, paste0("x", 1:3)]))
  sc <- (Z %*% p$chain$W[paste0("x", 1:3), "ind60"] - p$chain$score_center["ind60"]) / p$chain$score_scale["ind60"]
  expect_equal(unname(as.vector(sc)), unname(m$construct_scores[, "ind60"]), tolerance = 1e-8)
})

test_that("PLSc chain over-dispersion matches 1/sqrt(rho_A) on the single-antecedent mobi case", {
  # Tutorial oracle: on the full sample, SD(chain) / SD(implied) = 1.472 against 1/sqrt(.460) = 1.474.
  m <- quiet_pls(mobi, mobi_mm_plsc, mobi_sm)
  p <- sem_params(m)
  ch <- predict_oos(p, mobi, mobi_y, mobi_x, "chain")
  im <- predict_oos(p, mobi, mobi_y, mobi_x, "implied")
  ratio <- apply(ch, 2, sd) / apply(im, 2, sd)
  expect_equal(unname(ratio), rep(1.472, 3), tolerance = 2e-3)
  expect_equal(unname(1 / sqrt(p$rhoA["EXP"])), 1.474, tolerance = 2e-3)
})

test_that("earliest-antecedent chain propagates through the mediator", {
  m <- quiet_pls(pd, pd_mm_pls, pd_sm)
  p <- sem_params(m)
  ea <- predict_oos(p, pd, pd_y, pd_x_ea, "chain")
  expect_equal(dim(ea), c(75, 4))
  expect_false(anyNA(ea))
  # the EA prediction is a linear function of the ind60 score alone: rank one
  expect_equal(qr(scale(ea, scale = FALSE))$rank, 1)
})

test_that("predictions are on the raw scale (means close to the sample means)", {
  m <- quiet_pls(mobi, mobi_mm_pls, mobi_sm)
  p <- sem_params(m)
  im <- predict_oos(p, mobi, mobi_y, mobi_x, "implied")
  expect_equal(unname(colMeans(im)), unname(colMeans(mobi[, mobi_y])), tolerance = 1e-8)
})

test_that("unknown indicators and inadmissible params are refused", {
  m <- quiet_pls(mobi, mobi_mm_pls, mobi_sm)
  p <- sem_params(m)
  expect_error(predict_oos(p, mobi, "nope", mobi_x), "not in the fitted model")
  p$admissible <- FALSE; p$reason <- "test"
  expect_error(predict_oos(p, mobi, mobi_y, mobi_x, "implied"), "unavailable")
})
