#' Out-of-sample predictions from a fitted SEM
#'
#' Turns a fitted model into case-level predictions of observed outcome
#' indicators from observed predictor indicators, by one of two
#' constructions.
#'
#' * `"implied"`: the model-implied regression (de Rooij et al., 2023), the
#'   best linear prediction of `ynames` from `xnames` under the fitted mean
#'   vector and covariance matrix, `mu_y + Sigma_yx Sigma_xx^-1 (x - mu_x)`.
#'   For lavaan models this equals `lavaan::lavPredictY()`. For PLSc models
#'   the implied covariances are built from the corrected loadings and the
#'   corrected construct correlations, so scores and parameters are in the
#'   same metric; it is unavailable when the corrected solution is
#'   inadmissible (see [sem_params()]).
#' * `"chain"`: the construct-score chain used by PLSpredict (Shmueli et al.,
#'   2016, 2019): predictor items are weighted into construct scores, carried
#'   through the path coefficients in causal order, and mapped to the outcome
#'   items through the loadings. Only available for seminr models. Under
#'   consistent PLS this construction applies factor-metric coefficients to
#'   composite-metric scores and over-disperses the predictions.
#'
#' @param params An object from [sem_params()].
#' @param newdata Data frame containing `xnames`.
#' @param ynames,xnames Character vectors of outcome and predictor indicator
#'   names.
#' @param construction `"implied"` or `"chain"`.
#' @return A numeric matrix, rows of `newdata` by `ynames`.
#' @references de Rooij, M., et al. (2023). SEM-based out-of-sample
#'   predictions. *Structural Equation Modeling*, 30(1), 132–148.
#'   Shmueli, G., et al. (2016). The elephant in the room: Predictive
#'   performance of PLS models. *Journal of Business Research*, 69(10),
#'   4552–4564.
#' @examples
#' d <- lavaan::PoliticalDemocracy
#' fit <- lavaan::sem("ind60 =~ x1 + x2 + x3
#'                     dem60 =~ y1 + y2 + y3 + y4
#'                     dem65 =~ y5 + y6 + y7 + y8
#'                     dem60 ~ ind60
#'                     dem65 ~ ind60 + dem60", data = d[1:60, ], meanstructure = TRUE)
#' yhat <- predict_oos(sem_params(fit), d[61:75, ], ynames = paste0("y", 5:8),
#'                     xnames = c(paste0("x", 1:3), paste0("y", 1:4)))
#' head(yhat)
#' @export
predict_oos <- function(params, newdata, ynames, xnames,
                        construction = c("implied", "chain")) {
  stopifnot(inherits(params, "sem_params"))
  construction <- match.arg(construction)
  switch(construction,
         implied = predict_implied(params, newdata, ynames, xnames),
         chain = predict_chain(params, newdata, ynames, xnames))
}

predict_implied <- function(params, newdata, ynames, xnames) {
  if (!isTRUE(params$admissible))
    stop("model-implied regression unavailable: ", params$reason)
  S <- params$Sigma
  mu <- params$mu
  miss <- setdiff(c(ynames, xnames), colnames(S))
  if (length(miss))
    stop("not in the fitted model: ", paste(miss, collapse = ", "))
  Sxx <- S[xnames, xnames, drop = FALSE]
  Syx <- S[ynames, xnames, drop = FALSE]
  ev <- eigen(Sxx, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-10 * max(ev))
    stop("model-implied Sigma_xx is not positive definite")
  G <- Syx %*% solve(Sxx)
  X <- as.matrix(newdata[, xnames, drop = FALSE])
  Yhat <- sweep(X, 2, mu[xnames]) %*% t(G)
  Yhat <- sweep(Yhat, 2, mu[ynames], "+")
  dimnames(Yhat) <- list(NULL, ynames)
  Yhat
}

predict_chain <- function(params, newdata, ynames, xnames) {
  ch <- params$chain
  if (is.null(ch))
    stop("the construct-score chain is only available for seminr models")
  W <- ch$W; B <- ch$B; L <- ch$L; blocks <- ch$blocks; order <- ch$order
  miss <- setdiff(c(ynames, xnames), rownames(W))
  if (length(miss))
    stop("not in the fitted model: ", paste(miss, collapse = ", "))
  Z <- sweep(sweep(as.matrix(newdata[, xnames, drop = FALSE]), 2, params$mu[xnames]),
             2, params$sd[xnames], "/")
  eta <- matrix(NA_real_, nrow(Z), length(order), dimnames = list(NULL, order))
  for (cn in order) {
    items <- blocks[[cn]]
    if (all(items %in% xnames)) {
      sc <- Z[, items, drop = FALSE] %*% W[items, cn]
      eta[, cn] <- (sc - ch$score_center[cn]) / ch$score_scale[cn]
    } else {
      pred <- order[seq_len(match(cn, order) - 1)]
      pred <- pred[B[pred, cn] != 0]
      if (!length(pred)) next                     # exogenous and unobserved
      if (anyNA(eta[, pred])) next
      eta[, cn] <- eta[, pred, drop = FALSE] %*% B[pred, cn]
    }
  }
  Yz <- vapply(ynames, function(it) {
    cn <- names(blocks)[vapply(blocks, function(b) it %in% b, logical(1))][1]
    if (anyNA(eta[, cn]))
      stop("construct ", cn, " can be neither scored from `xnames` nor predicted from upstream constructs")
    eta[, cn] * L[it, cn]
  }, numeric(nrow(Z)))
  Yz <- matrix(Yz, ncol = length(ynames))
  Yhat <- sweep(sweep(Yz, 2, params$sd[ynames], "*"), 2, params$mu[ynames], "+")
  dimnames(Yhat) <- list(NULL, ynames)
  Yhat
}
