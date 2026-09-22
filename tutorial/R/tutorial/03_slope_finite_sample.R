# ==============================================================================
# 03_slope_finite_sample.R
# Finite-sample reference for the dispersion slope (supplement S3).
#
# Population: the ML-fitted PoliticalDemocracy model (no residual covariances),
# i.e. its implied mean vector and covariance matrix, so the population model is
# correct by construction and any departure of the slope from one is estimation
# error in theta-hat, not misfit.
#
# For each replication: draw n MVN cases from those moments; run 10-fold CV with
# lavaan + lavPredictY (the "estimated rule"), pool the out-of-sample predictions,
# and compute the centred observed-on-predicted slope; also compute the same slope
# for the ORACLE predictor that uses the true population moments (no estimation).
# Designs: DA (x1-x3, y1-y4 -> y5-y8) and EA (x1-x3 -> y5-y8).
# Run from the project root:  Rscript R/tutorial/03_slope_finite_sample.R
# ==============================================================================
suppressPackageStartupMessages({ library(lavaan); library(MASS) })
SEED <- 20260825; N_REP <- 200; K <- 10; NS <- c(75, 300, 1000)
out_dir <- "results/tutorial"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "03_slope_finite_sample.txt"), split = TRUE)
cat("Finite-sample dispersion-slope reference —", format(Sys.time()), "\n")
cat("lavaan", as.character(packageVersion("lavaan")), "| reps", N_REP, "| folds", K, "| seed", SEED, "\n\n")

lav_syntax <- '
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
  dem60 ~ ind60
  dem65 ~ ind60 + dem60
'
fit0 <- sem(lav_syntax, data = lavaan::PoliticalDemocracy, meanstructure = TRUE)
imp <- lavInspect(fit0, "implied")
Sigma <- imp$cov; mu <- imp$mean; vars <- rownames(Sigma)
y_names <- paste0("y", 5:8)
designs <- list(DA = c(paste0("x", 1:3), paste0("y", 1:4)), EA = paste0("x", 1:3))

slope_fn <- function(Y, Yhat) {           # pooled, centred within item
  yc <- scale(Y, scale = FALSE); hc <- scale(Yhat, scale = FALSE)
  sum(yc * hc) / sum(hc^2)
}
oracle_pred <- function(X, x_names) {
  G <- Sigma[y_names, x_names] %*% solve(Sigma[x_names, x_names])
  sweep(sweep(X, 2, mu[x_names]) %*% t(G), 2, mu[y_names], "+")
}

res <- list()
for (n in NS) {
  set.seed(SEED + n)
  est <- matrix(NA_real_, N_REP, length(designs), dimnames = list(NULL, names(designs)))
  ora <- est; nfail <- setNames(integer(length(designs)), names(designs))
  for (r in seq_len(N_REP)) {
    dat <- as.data.frame(mvrnorm(n, mu, Sigma)); names(dat) <- vars
    folds <- sample(rep(seq_len(K), length.out = n))
    fit_ok <- TRUE
    yhat <- setNames(vector("list", length(designs)), names(designs))
    for (d in names(designs)) yhat[[d]] <- matrix(NA_real_, n, length(y_names))
    for (k in seq_len(K)) {
      te <- folds == k
      f <- tryCatch(sem(lav_syntax, data = dat[!te, ], meanstructure = TRUE, warn = FALSE),
                    error = function(e) NULL)
      if (is.null(f) || !lavInspect(f, "converged")) { fit_ok <- FALSE; break }
      for (d in names(designs))
        yhat[[d]][te, ] <- as.matrix(lavPredictY(f, newdata = dat[te, ], ynames = y_names,
                                                xnames = designs[[d]]))
    }
    Y <- as.matrix(dat[, y_names])
    for (d in names(designs)) {
      if (fit_ok) est[r, d] <- slope_fn(Y, yhat[[d]]) else nfail[d] <- nfail[d] + 1L
      ora[r, d] <- slope_fn(Y, oracle_pred(as.matrix(dat[, designs[[d]]]), designs[[d]]))
    }
  }
  cat(sprintf("n = %4d  (replications with a non-converged fold, dropped: DA %d, EA %d)\n", n, nfail["DA"], nfail["EA"]))
  tab <- data.frame(design = names(designs),
                    estimated_rule = sprintf("%.3f (%.3f)", colMeans(est, na.rm = TRUE), apply(est, 2, sd, na.rm = TRUE)),
                    oracle_true_moments = sprintf("%.3f (%.3f)", colMeans(ora), apply(ora, 2, sd)))
  print(tab, row.names = FALSE); cat("\n")
  res[[as.character(n)]] <- list(est = est, ora = ora, nfail = nfail)
}
cat("Entries are mean (SD across replications) of the pooled centred slope of observed y5-y8 on predicted.\n")
pd <- file.path(out_dir, "01_political_democracy.rds")   # like-for-like: the reference above uses every fold
if (file.exists(pd)) { tb <- readRDS(pd)$tables
  cat(sprintf("Observed all-fold ML slopes (n = 75, supplement S6): DA %.3f, EA %.3f.\n",
              tb$DA$dispersion_slope[tb$DA$method == "cbsem"], tb$EA$dispersion_slope[tb$EA$method == "cbsem"])) }
sink()
saveRDS(res, file.path(out_dir, "03_slope_finite_sample.rds"))
