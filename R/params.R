#' Prediction parameters of a fitted SEM
#'
#' Extracts from a fitted model what the two prediction constructions need:
#' the model-implied indicator mean vector and covariance matrix (for the
#' model-implied regression) and, for `seminr` models, the outer weights, path
#' coefficients and loadings (for the construct-score chain). For consistent
#' PLS (constructs specified with `seminr::reflective()`) the construct
#' correlations are corrected for unreliability with each construct's
#' rho_A, propagated through the structural model, and checked for
#' admissibility: a corrected correlation of absolute value one or more, a
#' non-positive-definite construct correlation matrix, or a rho_A outside
#' (0, 1] makes the model-implied regression unavailable. The chain is
#' always available because it never forms the implied covariances.
#'
#' @param fit A fitted `lavaan` object (single group, fitted with
#'   `meanstructure = TRUE`) or a `seminr_model`.
#' @param ... Passed to methods.
#' @return An object of class `sem_params`: a list with `estimator`,
#'   `Sigma` and `mu` (model-implied moments on the raw indicator scale,
#'   `NULL` when inadmissible), `sd` (training standard deviations), `chain`
#'   (`NULL` for lavaan), `rhoA`, `admissible` (logical) and `reason`.
#' @examples
#' fit <- lavaan::sem("ind60 =~ x1 + x2 + x3
#'                     dem60 =~ y1 + y2 + y3 + y4
#'                     dem65 =~ y5 + y6 + y7 + y8
#'                     dem60 ~ ind60
#'                     dem65 ~ ind60 + dem60",
#'                    data = lavaan::PoliticalDemocracy, meanstructure = TRUE)
#' p <- sem_params(fit)
#' p$admissible
#' @export
sem_params <- function(fit, ...) UseMethod("sem_params")

#' @rdname sem_params
#' @export
sem_params.lavaan <- function(fit, ...) {
  if (!isTRUE(lavaan::lavInspect(fit, "converged")))
    stop("lavaan model did not converge")
  imp <- lavaan::lavInspect(fit, "implied")
  if (is.null(imp$cov))
    stop("only single-group lavaan models are supported")
  if (is.null(imp$mean))
    stop("fit the lavaan model with `meanstructure = TRUE`")
  Sigma <- unclass(as.matrix(imp$cov))
  mu <- as.numeric(imp$mean)
  names(mu) <- colnames(Sigma)
  new_sem_params(estimator = "lavaan", Sigma = Sigma, mu = mu,
                 sd = sqrt(diag(Sigma)), chain = NULL, rhoA = NULL,
                 admissible = TRUE, reason = NULL)
}

#' @rdname sem_params
#' @param plsc Logical; treat the model as consistent PLS. Defaults to
#'   `TRUE` when any construct was specified with `seminr::reflective()`.
#'   Composite constructs in a mixed model are given rho_A = 1.
#' @export
sem_params.seminr_model <- function(fit, plsc = NULL, ...) {
  W <- fit$outer_weights
  L <- fit$outer_loadings
  B <- fit$path_coef
  constructs <- colnames(W)
  items <- rownames(W)
  blocks <- lapply(stats::setNames(constructs, constructs),
                   function(cn) items[W[, cn] != 0])
  order <- topo_order(B)
  exo <- constructs[colSums(B != 0) == 0]
  types <- fit$mmMatrix[, "type"]
  reflective <- unique(fit$mmMatrix[types == "C", "construct"])
  if (is.null(plsc)) plsc <- length(reflective) > 0

  # training moments on the raw scale, and the training standardisation of
  # each composite score (seminr standardises scores on the training sample)
  X <- as.matrix(fit$data[, items, drop = FALSE])
  mu <- colMeans(X)
  s <- apply(X, 2, stats::sd)
  s[!is.finite(s) | s < 1e-10] <- 1
  Z <- sweep(sweep(X, 2, mu), 2, s, "/")
  sc <- vapply(constructs, function(cn) {
    v <- Z[, blocks[[cn]], drop = FALSE] %*% W[blocks[[cn]], cn]
    c(mean(v), stats::sd(v))
  }, numeric(2))
  chain <- list(W = W, B = B, L = L, blocks = blocks, order = order,
                score_center = sc[1, ], score_scale = sc[2, ])

  rA <- stats::setNames(rep(1, length(constructs)), constructs)
  Phi_s <- stats::cor(fit$construct_scores)[constructs, constructs]
  estimator <- if (plsc) "PLSc" else "PLS"
  bad <- function(msg)
    new_sem_params(estimator = estimator, Sigma = NULL, mu = mu, sd = s,
                   chain = chain, rhoA = rA, admissible = FALSE, reason = msg)
  if (plsc) {
    for (cn in reflective) rA[cn] <- as.numeric(seminr::rho_A(fit, cn))
    if (any(!is.finite(rA)) || any(rA <= 0) || any(rA > 1))
      return(bad(sprintf("inadmissible rho_A: %s",
                         paste(sprintf("%s = %.3f", constructs, rA), collapse = ", "))))
    Phi_s <- Phi_s / sqrt(outer(rA, rA))
    diag(Phi_s) <- 1
    off <- abs(Phi_s[upper.tri(Phi_s)])
    if (length(off) && max(off) >= 1)
      return(bad(sprintf("inadmissible: disattenuated construct correlation >= 1 (max %.3f)",
                         max(off))))
  }
  Phi <- tryCatch(implied_construct_cor(B, Phi_s, exo, order), error = function(e) e)
  if (inherits(Phi, "error")) return(bad(conditionMessage(Phi)))
  if (min(eigen(Phi, symmetric = TRUE, only.values = TRUE)$values) <= 1e-8)
    return(bad("implied construct correlation matrix not positive definite"))
  Theta <- diag(pmax(1 - rowSums(L^2), 1e-6), nrow = length(items))
  dimnames(Theta) <- list(items, items)
  Sigma_z <- L %*% Phi %*% t(L) + Theta
  Sigma <- Sigma_z * outer(s, s)
  dimnames(Sigma) <- list(items, items)
  out <- new_sem_params(estimator = estimator, Sigma = Sigma, mu = mu, sd = s,
                        chain = chain, rhoA = rA, admissible = TRUE, reason = NULL)
  out$Lambda <- L
  out$Phi <- Phi
  out$Theta <- Theta
  out
}

new_sem_params <- function(estimator, Sigma, mu, sd, chain, rhoA, admissible, reason) {
  structure(list(estimator = estimator, Sigma = Sigma, mu = mu, sd = sd,
                 chain = chain, rhoA = rhoA, admissible = admissible,
                 reason = reason), class = "sem_params")
}

#' @export
print.sem_params <- function(x, ...) {
  cat("<sem_params>", x$estimator, "\n")
  cat("  indicators:", length(x$mu), "\n")
  if (!is.null(x$chain))
    cat("  constructs (causal order):", paste(x$chain$order, collapse = " -> "), "\n")
  if (!is.null(x$rhoA))
    cat("  rho_A:", paste(sprintf("%s %.3f", names(x$rhoA), x$rhoA), collapse = "; "), "\n")
  cat("  model-implied regression available:", x$admissible,
      if (!x$admissible) paste0(" (", x$reason, ")") else "", "\n")
  invisible(x)
}

# Structural-model-implied construct correlation matrix for a recursive model.
# B[k, j] is the path from k to j (seminr layout). Phi_exo holds the
# correlations among exogenous constructs. Constructs are standardised, so the
# disturbance variance of each endogenous construct is whatever makes its
# variance one.
implied_construct_cor <- function(B, Phi_exo, exo, order) {
  p <- length(order)
  Phi <- matrix(0, p, p, dimnames = list(order, order))
  Phi[exo, exo] <- Phi_exo[exo, exo]
  for (j in setdiff(order, exo)) {
    pred <- order[seq_len(match(j, order) - 1)]
    b <- B[pred, j]
    cov_j_pred <- as.numeric(Phi[pred, pred, drop = FALSE] %*% b)
    Phi[j, pred] <- Phi[pred, j] <- cov_j_pred
    Phi[j, j] <- 1
  }
  Phi
}

# Causal (topological) order of the constructs from the path matrix.
topo_order <- function(B) {
  constructs <- colnames(B)
  remaining <- constructs
  out <- character(0)
  while (length(remaining)) {
    free <- remaining[colSums(B[remaining, remaining, drop = FALSE] != 0) == 0]
    if (!length(free)) stop("structural model is not recursive")
    out <- c(out, free)
    remaining <- setdiff(remaining, free)
  }
  out
}
