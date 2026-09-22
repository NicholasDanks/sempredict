# ==============================================================================
# 00_helpers.R — shared machinery for the tutorial scripts
#
# General (any recursive model) implementations of:
#   implied_construct_cor()  structural-model-implied construct correlations
#   predict_chain()          composite-score chain  W -> B -> Lambda  (PLSpredict)
#   predict_cov_blp()       covariance-based best linear prediction (de Rooij 2023)
#   seminr_params()          (W, B, Lambda, Phi, Theta) from a seminr fit, with the
#                            PLSc admissibility check
#   run_design() / summarise_design()   repeated k-fold CV, item RMSE, dispersion
#                            slope, OFR, descriptive paired-loss differences, failure counts
#
# The CV functions read these objects from the calling environment (deliberately,
# to keep the tutorial scripts linear): data, n, designs, blocks, causal_order,
# exo_constructs, mm_pls, mm_plsc, sm, lav_syntax, methods, N_REPS, K_FOLDS, SEED.
# ==============================================================================

# ---- General helpers ---------------------------------------------------------

standardise <- function(newdata, train_data, vars) {
  m <- colMeans(train_data[, vars, drop = FALSE])
  s <- apply(train_data[, vars, drop = FALSE], 2, sd); s[s < 1e-10] <- 1
  list(z = sweep(sweep(as.matrix(newdata[, vars, drop = FALSE]), 2, m, "-"), 2, s, "/"),
       m = m, s = s)
}

# Structural-model-implied construct correlation matrix for a recursive model.
# B[k, j] = path from k to j (seminr path_coef layout). Phi_exo = correlations
# among exogenous constructs. Disturbance variances chosen so diag(Phi) = 1
# (standardised constructs), processed in causal order.
implied_construct_cor <- function(B, Phi_exo, exo, order) {
  p <- length(order); Phi <- matrix(0, p, p, dimnames = list(order, order))
  Phi[exo, exo] <- Phi_exo[exo, exo]
  for (j in setdiff(order, exo)) {
    pred <- order[seq_len(match(j, order) - 1)]
    b <- B[pred, j]
    cov_j_pred <- as.numeric(Phi[pred, pred, drop = FALSE] %*% b)   # Cov(eta_j, pred)
    Phi[j, pred] <- Phi[pred, j] <- cov_j_pred
    Phi[j, j] <- 1   # disturbance variance = 1 - b' Phi_pred b (standardised)
  }
  Phi
}

# Composite-score chain (PLSpredict logic), general for recursive models.
# Constructs whose indicators are all in x_names get scores from data (W);
# every other construct is propagated through B in causal order. With
# x = direct antecedents this is predict_DA; with x = earliest antecedent only
# it is predict_EA. Loadings map the endogenous score to its indicators.
predict_chain <- function(W, B, L, newdata, train_data, x_names, y_names,
                          blocks, order) {
  z <- standardise(newdata, train_data, x_names)$z
  eta <- matrix(NA_real_, nrow(newdata), length(order), dimnames = list(NULL, order))
  for (cn in order) {
    items <- blocks[[cn]]
    if (all(items %in% x_names)) {
      w <- W[items, cn]
      sc <- z[, items, drop = FALSE] %*% w
      # seminr standardises scores on the training sample
      tr <- standardise(train_data, train_data, items)$z %*% w
      eta[, cn] <- (sc - mean(tr)) / sd(tr)
    } else {
      pred <- order[seq_len(match(cn, order) - 1)]
      eta[, cn] <- eta[, pred, drop = FALSE] %*% B[pred, cn]
    }
  }
  st <- standardise(newdata, train_data, y_names)
  Yz <- sapply(y_names, function(it) {
    cn <- names(blocks)[sapply(blocks, function(b) it %in% b)]
    eta[, cn] * L[it, cn]
  })
  sweep(sweep(matrix(Yz, ncol = length(y_names)), 2, st$s, "*"), 2, st$m, "+") |>
    `colnames<-`(y_names)
}

# Covariance-based best linear prediction (de Rooij et al. 2023) from
# (Lambda, Phi, Theta) in the standardised metric, mapped back to raw scale.
predict_cov_blp <- function(L, Phi, Theta, newdata, train_data, x_names, y_names) {
  Sigma <- L %*% Phi %*% t(L) + Theta
  Sxx <- Sigma[x_names, x_names, drop = FALSE]
  Syx <- Sigma[y_names, x_names, drop = FALSE]
  ev <- eigen(Sxx, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-8) stop("Sigma_xx not positive definite (inadmissible)")
  Gamma <- Syx %*% solve(Sxx)
  z <- standardise(newdata, train_data, x_names)$z
  st <- standardise(newdata, train_data, y_names)
  Yz <- z %*% t(Gamma)
  sweep(sweep(Yz, 2, st$s, "*"), 2, st$m, "+") |> `colnames<-`(y_names)
}

# Extract (W, B, L, Phi, Theta) from a seminr fit. For PLSc (reflective()),
# B and L are already in the factor metric; Phi is the score correlation
# disattenuated by rho_A, then constrained by the structural model. For plain
# PLS the score correlations are used as they are.
# Admissibility is a property of Sigma(theta), which only the covariance-based prediction
# uses. W, B and L are always returned so that the composite-score chain -- which
# never forms Phi, Theta or rho_A -- can be evaluated on every fold. An
# inadmissible fit is reported in $inadmissible rather than raised, so the caller
# can fail the rule alone and keep the two rows comparable on the folds they share.
seminr_params <- function(m, plsc, blocks, order, exo) {
  W <- m$outer_weights; L <- m$outer_loadings; B <- m$path_coef
  rA <- rep(NA_real_, length(order)); names(rA) <- order
  l2 <- rowSums(L^2)                        # > 1 means a negative residual variance
  also <- if (plsc && any(l2 > 1)) "; a loading also exceeds one" else ""
  bad <- function(msg) list(W = W, B = B, L = L, Phi = NULL, Theta = NULL,
                            rhoA = rA, inadmissible = paste0(msg, also))
  Phi_s <- cor(m$construct_scores)[order, order]
  if (plsc) {
    rA <- sapply(order, function(cn) seminr::rho_A(m, cn))
    if (any(!is.finite(rA)) || any(rA <= 0) || any(rA > 1))
      return(bad(sprintf("inadmissible rho_A: %s", paste(round(rA, 3), collapse = ", "))))
    Phi_s <- Phi_s / sqrt(outer(rA, rA))
    diag(Phi_s) <- 1
    if (max(abs(Phi_s[upper.tri(Phi_s)])) >= 1)
      return(bad(sprintf("inadmissible: disattenuated construct correlation >= 1 (max %.3f)",
                         max(abs(Phi_s[upper.tri(Phi_s)])))))
  }
  Phi <- tryCatch(implied_construct_cor(B, Phi_s, exo, order), error = function(e) e)
  if (inherits(Phi, "error")) return(bad(conditionMessage(Phi)))
  if (min(eigen(Phi, symmetric = TRUE, only.values = TRUE)$values) <= 1e-8)
    return(bad("implied construct correlation matrix not positive definite"))
  # A standardised loading above one implies a negative residual variance (a Heywood case):
  # no covariance matrix of real data has one, so the solution is inadmissible. It is NOT
  # repaired by flooring the residual variance: a floored item becomes an error-free measure
  # of its construct and takes nearly all the weight in Sigma_xx^{-1}.
  if (any(l2 > 1)) { also <- ""
    return(bad(sprintf("inadmissible: standardised loading above one (%s)",
                       paste(names(l2)[l2 > 1], collapse = ", ")))) }
  Theta <- diag(1 - l2); dimnames(Theta) <- dimnames(L)[c(1, 1)]
  list(W = W, B = B, L = L, Phi = Phi, Theta = Theta, rhoA = rA, inadmissible = NULL)
}

# Linear prediction map of each method on standardised x (Gamma: y-items x x-items),
# so the credit given to each antecedent block can be compared directly.
gamma_chain <- function(par, x_names, y_names, blocks, order) {
  # eta_exo = W' z (then standardised on training scores -- ignore the scale for
  # the comparison: seminr scores have unit variance, W already normalised)
  p <- length(x_names); G <- matrix(0, length(y_names), p, dimnames = list(y_names, x_names))
  eta_map <- matrix(0, length(order), p, dimnames = list(order, x_names))
  for (cn in order) {
    items <- blocks[[cn]]
    if (all(items %in% x_names)) eta_map[cn, items] <- par$W[items, cn]
    else { pred <- order[seq_len(match(cn, order) - 1)]
           eta_map[cn, ] <- as.numeric(t(par$B[pred, cn, drop = FALSE]) %*% eta_map[pred, , drop = FALSE]) }
  }
  for (it in y_names) { cn <- names(blocks)[sapply(blocks, function(b) it %in% b)]
                        G[it, ] <- par$L[it, cn] * eta_map[cn, ] }
  G
}
gamma_sem <- function(par, x_names, y_names) {
  Sigma <- par$L %*% par$Phi %*% t(par$L) + par$Theta
  Sigma[y_names, x_names] %*% solve(Sigma[x_names, x_names])
}
block_credit <- function(G, blocks, x_names) {
  # summed absolute prediction weight per antecedent block, averaged over y items
  sapply(names(blocks)[sapply(blocks, function(b) all(b %in% x_names))],
         function(cn) mean(rowSums(abs(G[, blocks[[cn]], drop = FALSE]))))
}

# TRUE when a params object supports the covariance-based prediction (W/B/L alone are
# enough for the chain, so the chain must not be gated on this).
params_ok <- function(p) !is.null(p) && !inherits(p, "error") && is.null(p$inadmissible)

fit_lavaan <- function(train) {
  fit <- tryCatch(sem(lav_syntax, data = train, meanstructure = TRUE,
                      warn = FALSE), error = function(e) NULL)
  if (is.null(fit) || !lavInspect(fit, "converged")) return(NULL)
  fit
}

# ---- Section B: repeated k-fold CV ------------------------------------------
methods <- c("cbsem", "pls_chain", "pls_cov_plugin", "plsc_chain", "plsc_cov_blp", "lm", "mean")

run_design <- function(dname) {
  x_names <- designs[[dname]]$x; y_names <- designs[[dname]]$y
  # per-observation, per-rep, per-method squared error (averaged over y items)
  se_obs <- array(NA_real_, c(n, N_REPS, length(methods)), dimnames = list(NULL, NULL, methods))
  se_item <- array(NA_real_, c(n, N_REPS, length(methods), length(y_names)),
                   dimnames = list(NULL, NULL, methods, y_names))
  yhat <- array(NA_real_, c(n, N_REPS, length(methods), length(y_names)),
                dimnames = list(NULL, NULL, methods, y_names))
  mse_in <- array(NA_real_, c(N_REPS, K_FOLDS, length(methods)), dimnames = list(NULL, NULL, methods))
  fail <- setNames(integer(length(methods)), methods)
  fail_reason <- character(0)
  rhoA_log <- list()

  set.seed(SEED)
  for (r in seq_len(N_REPS)) {
    folds <- sample(rep(seq_len(K_FOLDS), length.out = n))
    for (k in seq_len(K_FOLDS)) {
      te <- folds == k; train <- data[!te, ]; test <- data[te, ]
      Ytr <- as.matrix(train[, y_names]); Yte <- as.matrix(test[, y_names])

      fits <- list()
      fits$lav  <- fit_lavaan(train)
      fits$pls  <- tryCatch(estimate_pls(train, mm_pls,  sm), error = function(e) NULL)
      fits$plsc <- tryCatch(estimate_pls(train, mm_plsc, sm), error = function(e) NULL)
      par <- list()
      par$pls  <- if (!is.null(fits$pls))  tryCatch(seminr_params(fits$pls,  FALSE, blocks, causal_order, exo_constructs), error = function(e) e) else NULL
      par$plsc <- if (!is.null(fits$plsc)) tryCatch(seminr_params(fits$plsc, TRUE,  blocks, causal_order, exo_constructs), error = function(e) e) else NULL
      if (!is.null(fits$plsc)) rhoA_log[[length(rhoA_log) + 1]] <-
        sapply(causal_order, function(cn) seminr::rho_A(fits$plsc, cn))

      pred <- function(method, newdata) {
        switch(method,
          cbsem = if (is.null(fits$lav)) stop("lavaan did not converge") else
                    lavPredictY(fits$lav, newdata = newdata, ynames = y_names, xnames = x_names),
          pls_chain  = with(par$pls,  predict_chain(W, B, L, newdata, train, x_names, y_names, blocks, causal_order)),
          pls_cov_plugin    = with(par$pls,  predict_cov_blp(L, Phi, Theta, newdata, train, x_names, y_names)),
          plsc_chain = with(par$plsc, predict_chain(W, B, L, newdata, train, x_names, y_names, blocks, causal_order)),
          plsc_cov_blp   = with(par$plsc, predict_cov_blp(L, Phi, Theta, newdata, train, x_names, y_names)),
          lm   = { f <- lm(as.formula(sprintf("cbind(%s) ~ %s", paste(y_names, collapse = ","),
                                               paste(x_names, collapse = "+"))), data = train)
                   predict(f, newdata = newdata) },
          mean = matrix(colMeans(Ytr), nrow(newdata), length(y_names), byrow = TRUE,
                        dimnames = list(NULL, y_names)))
      }
      for (mth in methods) {
        if (grepl("^plsc", mth)) {
          # PLSc estimation itself failed: neither rule nor chain is defined.
          if (is.null(par$plsc) || inherits(par$plsc, "error")) { fail[mth] <- fail[mth] + 1L
            fail_reason <- c(fail_reason, if (is.null(par$plsc)) "PLSc estimation did not converge"
                                          else conditionMessage(par$plsc)); next }
          # Sigma(theta) inadmissible: the covariance-based prediction is undefined, the chain is not.
          if (mth == "plsc_cov_blp" && !is.null(par$plsc$inadmissible)) { fail[mth] <- fail[mth] + 1L
            fail_reason <- c(fail_reason, par$plsc$inadmissible); next }
        }
        out <- tryCatch(pred(mth, test), error = function(e) e)
        if (inherits(out, "error")) { fail[mth] <- fail[mth] + 1L; fail_reason <- c(fail_reason, conditionMessage(out)); next }
        out <- as.matrix(out)[, y_names, drop = FALSE]
        se_item[te, r, mth, ] <- (Yte - out)^2
        yhat[te, r, mth, ] <- out
        se_obs[te, r, mth] <- rowMeans((Yte - out)^2)
        ins <- tryCatch(as.matrix(pred(mth, train))[, y_names, drop = FALSE], error = function(e) NULL)
        if (!is.null(ins)) mse_in[r, k, mth] <- mean((Ytr - ins)^2)
      }
    }
  }
  list(design = dname, x = x_names, y = y_names, se_obs = se_obs, se_item = se_item, yhat = yhat,
       mse_in = mse_in, fail = fail, fail_reason = fail_reason,
       rhoA = do.call(rbind, rhoA_log))
}

summarise_design <- function(res) {
  cat(sprintf("\n--- Design %s: predict %s from %s ---\n", res$design,
              paste(res$y, collapse = ","), paste(res$x, collapse = ",")))
  keep <- methods
  # Per-repetition item RMSE, averaged across outcome indicators.
  rep_rmse <- sapply(keep, function(m) sapply(seq_len(N_REPS), function(r) {
    mean(sapply(seq_along(res$y), function(j)
      sqrt(mean(res$se_item[, r, m, j], na.rm = TRUE))), na.rm = TRUE)
  }))
  rmse <- colMeans(rep_rmse, na.rm = TRUE)
  rmse_mcse <- apply(rep_rmse, 2, sd, na.rm = TRUE) / sqrt(N_REPS)

  mse_out <- apply(res$se_obs, 3, mean, na.rm = TRUE)
  mse_in <- apply(res$mse_in, 3, mean, na.rm = TRUE)
  ofr <- (mse_out - mse_in) / mse_in

  # Descriptive paired loss differences only. Repeated k-fold training sets overlap,
  # so the tutorial does not attach t-test p values to these repeated-CV summaries.
  loss <- apply(res$se_obs, c(1, 3), mean, na.rm = TRUE)
  loss_diff <- function(m, bench) mean(loss[, bench] - loss[, m], na.rm = TRUE)

  # Pooled centered dispersion slope, computed within repetition then averaged.
  Yobs <- as.matrix(data[, res$y])
  rep_slope <- sapply(keep, function(m) sapply(seq_len(N_REPS), function(r) {
    yh <- res$yhat[, r, m, ]; if (all(is.na(yh))) return(NA_real_)
    ok <- complete.cases(yh)
    yc <- scale(Yobs[ok, , drop = FALSE], scale = FALSE)
    hc <- scale(yh[ok, , drop = FALSE], scale = FALSE)
    den <- sum(hc^2)
    if (!is.finite(den) || den <= 1e-12) return(NA_real_)
    sum(yc * hc) / den
  }))
  slope <- colMeans(rep_slope, na.rm = TRUE)
  slope_mcse <- apply(rep_slope, 2, sd, na.rm = TRUE) / sqrt(N_REPS)

  tab <- data.frame(
    method = keep,
    # two more decimals than the article reports, so article values are rounded once, not twice
    RMSE = round(rmse[keep], 5), RMSE_MCSE = round(rmse_mcse[keep], 5),
    dispersion_slope = round(slope[keep], 4), slope_MCSE = round(slope_mcse[keep], 4),
    OFR = round(ofr[keep], 4), failed_folds = res$fail[keep],
    dLoss_vs_LM = round(sapply(keep, loss_diff, bench = "lm"), 3),
    dLoss_vs_mean = round(sapply(keep, loss_diff, bench = "mean"), 3)
  )
  print(tab, row.names = FALSE)
  if (any(res$fail > 0)) {
    cat("\nFailure reasons (count):\n"); print(table(res$fail_reason))
  }
  if (!is.null(res$rhoA)) {
    cat(sprintf("\nrho_A across %d PLSc training fits: ", nrow(res$rhoA)))
    cat(paste(sprintf("%s [%.2f, %.2f]", colnames(res$rhoA),
                      apply(res$rhoA, 2, min), apply(res$rhoA, 2, max)), collapse = "; "), "\n")
  }
  invisible(tab)
}

# Rescore every method on the folds for which the PLSc covariance BLP is available.
# This produces the like-for-like comparison reported in Table 2 of the article.
summarise_common_plsc_folds <- function(res) {
  keep <- c("cbsem", "pls_chain", "plsc_chain", "plsc_cov_blp", "lm", "mean")
  common <- is.finite(res$se_obs[, , "plsc_cov_blp"])
  Yobs <- as.matrix(data[, res$y])

  # Recreate the fold assignments used by run_design() so in-sample MSE is
  # restricted to the same admissible training/test folds.
  set.seed(SEED)
  fold_id <- matrix(NA_integer_, n, N_REPS)
  for (r in seq_len(N_REPS))
    fold_id[, r] <- sample(rep(seq_len(K_FOLDS), length.out = n))

  rep_rmse <- matrix(NA_real_, N_REPS, length(keep), dimnames = list(NULL, keep))
  rep_slope <- rep_rmse
  rep_ofr <- rep_rmse

  for (r in seq_len(N_REPS)) {
    ok <- common[, r]
    ok_folds <- sapply(seq_len(K_FOLDS), function(k) all(ok[fold_id[, r] == k]))
    for (m in keep) {
      rep_rmse[r, m] <- mean(sapply(seq_along(res$y), function(j)
        sqrt(mean(res$se_item[ok, r, m, j], na.rm = TRUE))), na.rm = TRUE)

      yh <- res$yhat[ok, r, m, , drop = FALSE]
      yh <- matrix(yh, nrow = sum(ok), ncol = length(res$y))
      if (m != "mean") {
        yc <- scale(Yobs[ok, , drop = FALSE], scale = FALSE)
        hc <- scale(yh, scale = FALSE)
        den <- sum(hc^2)
        if (is.finite(den) && den > 1e-12) rep_slope[r, m] <- sum(yc * hc) / den
      }

      mse_out_r <- mean(res$se_obs[ok, r, m], na.rm = TRUE)
      mse_in_r <- mean(res$mse_in[r, ok_folds, m], na.rm = TRUE)
      rep_ofr[r, m] <- (mse_out_r - mse_in_r) / mse_in_r
    }
  }

  data.frame(
    method = keep,
    RMSE = colMeans(rep_rmse, na.rm = TRUE),
    RMSE_MCSE = apply(rep_rmse, 2, sd, na.rm = TRUE) / sqrt(N_REPS),
    dispersion_slope = colMeans(rep_slope, na.rm = TRUE),
    slope_MCSE = apply(rep_slope, 2, sd, na.rm = TRUE) / sqrt(N_REPS),
    OFR = colMeans(rep_ofr, na.rm = TRUE),
    common_obs_rep_cells = sum(common)
  )
}

