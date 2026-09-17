#' Summarise a cross-validated prediction run
#'
#' For every result column of a [cv_predict()] run: RMSE averaged over outcome
#' indicators, its Monte Carlo standard error across repetitions, the pooled
#' centred dispersion slope of observed on predicted (1 = ideal; below 1 the
#' predictions are over-dispersed), the overfit ratio
#' `OFR = (MSE_out - MSE_in) / MSE_in`, and the number of failed folds.
#'
#' The Monte Carlo standard error describes how much a number would move if
#' the folds were drawn again. It says nothing about variation across
#' samples, which is usually larger.
#'
#' @param cv An object from [cv_predict()].
#' @param common Optional name of a result column. When given, every column
#'   is rescored on the folds for which that column produced a prediction, so
#'   that a construction that is unavailable on some folds (an inadmissible
#'   PLSc solution, say) is compared like for like.
#' @param digits Rounding for the printed table; `NULL` for none.
#' @return A data frame with one row per result column and attribute
#'   `n_cells` (case-by-repetition cells scored).
#' @export
assess <- function(cv, common = NULL, digits = 4) {
  stopifnot(inherits(cv, "sem_cv"))
  cols <- cv$columns; reps <- cv$reps; k <- cv$k
  Yobs <- as.matrix(cv$data[, cv$ynames, drop = FALSE])
  ok_all <- matrix(TRUE, nrow(cv$data), reps)
  if (!is.null(common)) {
    if (!common %in% cols) stop("`common` must be one of: ", paste(cols, collapse = ", "))
    ok_all <- is.finite(cv$se_obs[, , common, drop = TRUE])
    ok_all <- matrix(ok_all, nrow(cv$data), reps)
  }
  rep_rmse <- rep_slope <- rep_ofr <- matrix(NA_real_, reps, length(cols), dimnames = list(NULL, cols))
  for (r in seq_len(reps)) {
    ok <- ok_all[, r]
    ok_folds <- vapply(seq_len(k), function(f) all(ok[cv$folds[, r] == f]), logical(1))
    for (m in cols) {
      se <- cv$se_item[ok, r, m, , drop = FALSE]
      rep_rmse[r, m] <- mean(vapply(seq_along(cv$ynames), function(j)
        sqrt(mean(se[, , , j], na.rm = TRUE)), numeric(1)), na.rm = TRUE)
      yh <- matrix(cv$yhat[ok, r, m, , drop = FALSE], nrow = sum(ok))
      keep <- stats::complete.cases(yh)
      if (any(keep) && m != "mean") {
        yc <- scale(Yobs[ok, , drop = FALSE][keep, , drop = FALSE], scale = FALSE)
        hc <- scale(yh[keep, , drop = FALSE], scale = FALSE)
        den <- sum(hc^2)
        if (is.finite(den) && den > 1e-12) rep_slope[r, m] <- sum(yc * hc) / den
      }
      mse_out <- mean(cv$se_obs[ok, r, m], na.rm = TRUE)
      mse_in <- mean(cv$mse_in[r, ok_folds, m], na.rm = TRUE)
      rep_ofr[r, m] <- (mse_out - mse_in) / mse_in
    }
  }
  mcse <- function(M) if (reps > 1) apply(M, 2, stats::sd, na.rm = TRUE) / sqrt(reps) else rep(NA_real_, ncol(M))
  out <- data.frame(
    method = cols,
    RMSE = colMeans(rep_rmse, na.rm = TRUE),
    RMSE_MCSE = mcse(rep_rmse),
    dispersion_slope = colMeans(rep_slope, na.rm = TRUE),
    slope_MCSE = mcse(rep_slope),
    OFR = colMeans(rep_ofr, na.rm = TRUE),
    failed_folds = unname(cv$fail[cols]),
    row.names = NULL, stringsAsFactors = FALSE)
  if (!is.null(digits)) for (v in c("RMSE", "RMSE_MCSE", "dispersion_slope", "slope_MCSE", "OFR"))
    out[[v]] <- round(out[[v]], digits)
  attr(out, "n_cells") <- sum(ok_all)
  attr(out, "common") <- common
  class(out) <- c("sem_assessment", "data.frame")
  out
}

#' @export
print.sem_assessment <- function(x, ...) {
  if (!is.null(attr(x, "common")))
    cat(sprintf("Rescored on the folds where '%s' was available (%d case-by-repetition cells)\n",
                attr(x, "common"), attr(x, "n_cells")))
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}

#' Failure log of a cross-validated run
#'
#' @param cv An object from [cv_predict()].
#' @return A data frame of result column, reason and count.
#' @export
failures <- function(cv) {
  stopifnot(inherits(cv, "sem_cv"))
  rows <- lapply(cv$columns, function(m) {
    if (!length(cv$reasons[[m]])) return(NULL)
    tb <- table(cv$reasons[[m]])
    data.frame(method = m, reason = names(tb), count = as.integer(tb), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) data.frame(method = character(0), reason = character(0), count = integer(0)) else out
}

#' Cross-validated predictive ability test (CVPAT)
#'
#' The paired test of Liengaard et al. (2021), applied as published on one
#' cross-validation repetition: each held-out case's squared error is
#' averaged over the outcome indicators, and a paired t-test compares the
#' benchmark's loss with the model's over the cases for which both are
#' finite. A positive difference favours the model.
#'
#' @param cv An object from [cv_predict()].
#' @param model,benchmark Result column names.
#' @param rep Which repetition to use (default the first).
#' @param alternative Passed to [stats::t.test()].
#' @return A one-row data frame: `model`, `benchmark`, `rep`, `n`, `d` (mean
#'   loss difference, benchmark minus model), `t`, `p`.
#' @references Liengaard, B. D., et al. (2021). Prediction: Coveted, yet
#'   forsaken? Introducing a cross-validated predictive ability test in
#'   partial least squares path modeling. *Decision Sciences*, 52(2), 362–392.
#' @export
cvpat <- function(cv, model, benchmark = "lm", rep = 1, alternative = "two.sided") {
  stopifnot(inherits(cv, "sem_cv"))
  for (v in c(model, benchmark)) if (!v %in% cv$columns) stop("unknown result column: ", v)
  d <- cv$se_obs[, rep, benchmark] - cv$se_obs[, rep, model]
  d <- d[is.finite(d)]
  if (length(d) < 3) return(data.frame(model = model, benchmark = benchmark, rep = rep,
                                       n = length(d), d = NA_real_, t = NA_real_, p = NA_real_))
  tt <- stats::t.test(d, alternative = alternative)
  data.frame(model = model, benchmark = benchmark, rep = rep, n = length(d),
             d = mean(d), t = unname(tt$statistic), p = tt$p.value, stringsAsFactors = FALSE)
}

#' CVPAT for every result column against one benchmark
#'
#' @inheritParams cvpat
#' @param benchmark Result column used as the benchmark.
#' @param models Result columns to test; defaults to all except the benchmark.
#' @param sensitivity If `TRUE` and the run has several repetitions, adds the
#'   range of `p` across repetitions. Repetitions share training data, so this
#'   is a sensitivity range, not independent evidence.
#' @return A data frame, one row per model.
#' @export
cvpat_table <- function(cv, benchmark = "lm", models = NULL, rep = 1, sensitivity = TRUE) {
  if (is.null(models)) models <- setdiff(cv$columns, benchmark)
  out <- do.call(rbind, lapply(models, function(m) cvpat(cv, m, benchmark, rep)))
  if (sensitivity && cv$reps > 1) {
    ps <- sapply(models, function(m) sapply(seq_len(cv$reps), function(r) cvpat(cv, m, benchmark, r)$p))
    ps <- matrix(ps, ncol = length(models))
    out$p_min <- apply(ps, 2, min, na.rm = TRUE)
    out$p_max <- apply(ps, 2, max, na.rm = TRUE)
  }
  out
}
