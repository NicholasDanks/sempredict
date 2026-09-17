#' Define a method for cross-validated prediction
#'
#' A method is a fitting function plus one or more prediction constructions.
#' [cv_predict()] refits it on every training fold and predicts the held-out
#' fold by each construction.
#'
#' @param fit A function of one argument, the training data frame, returning
#'   a fitted `lavaan` or `seminr_model` object. Errors and non-convergence
#'   are caught and counted as failures.
#' @param construction Character vector, any of `"implied"` and `"chain"`.
#'   With two constructions the method contributes two result columns,
#'   `label_implied` and `label_chain`, from one fit per fold.
#' @param label Name used in results; defaults to the name given in the
#'   `methods` list of [cv_predict()].
#' @param params_args List of extra arguments to [sem_params()] (for example
#'   `list(plsc = TRUE)`).
#' @return An object of class `sem_method`.
#' @examples
#' m <- sem_method(function(d) lavaan::sem("f =~ x1 + x2 + x3", data = d,
#'                                          meanstructure = TRUE))
#' @export
sem_method <- function(fit, construction = "implied", label = NULL, params_args = list()) {
  stopifnot(is.function(fit))
  construction <- match.arg(construction, c("implied", "chain"), several.ok = TRUE)
  structure(list(fit = fit, construction = construction, label = label,
                 params_args = params_args), class = "sem_method")
}

#' Repeated k-fold cross-validated prediction with identical folds
#'
#' Splits the data into `k` folds, `reps` times, and for every fold refits
#' every method on the training part and predicts `ynames` for the held-out
#' part from `xnames`. All methods and benchmarks see the same folds, so every
#' comparison is paired. Failures (non-convergence, inadmissible PLSc
#' solutions, prediction errors) are counted per method with their reasons
#' and never stop the run.
#'
#' @param methods Named list of [sem_method()] objects.
#' @param data Data frame with all indicators.
#' @param ynames,xnames Outcome and predictor indicator names.
#' @param k Number of folds.
#' @param reps Number of repetitions (fresh random fold assignments).
#' @param seed Random seed for the fold assignments; set it so the run is
#'   reproducible.
#' @param benchmarks Any of `"mean"` (the training-mean predictor) and `"lm"`
#'   (a multivariate linear regression of `ynames` on `xnames`).
#' @param verbose Print progress.
#' @return An object of class `sem_cv` holding the fold assignments, the
#'   out-of-sample predictions and squared errors for every case, repetition
#'   and method, the in-sample MSE per training fold, and the failure log.
#'   Summarise it with [assess()], test it with [cvpat()], draw it with
#'   [plot_dispersion()].
#' @examples
#' d <- seminr::mobi
#' mm <- seminr::constructs(
#'   seminr::reflective("EXP", seminr::multi_items("CUEX", 1:3)),
#'   seminr::reflective("SAT", seminr::multi_items("CUSA", 1:3)))
#' sm <- seminr::relationships(seminr::paths(from = "EXP", to = "SAT"))
#' plsc <- sem_method(function(tr) seminr::estimate_pls(tr, mm, sm),
#'                    construction = c("chain", "implied"))
#' cv <- cv_predict(list(PLSc = plsc), d, ynames = paste0("CUSA", 1:3),
#'                  xnames = paste0("CUEX", 1:3), k = 10, reps = 2, seed = 1)
#' assess(cv)
#' @export
cv_predict <- function(methods, data, ynames, xnames, k = 10, reps = 1, seed = NULL,
                       benchmarks = c("mean", "lm"), verbose = FALSE) {
  stopifnot(is.list(methods), length(methods) > 0, !is.null(names(methods)),
            all(nzchar(names(methods))))
  for (m in methods) if (!inherits(m, "sem_method")) stop("`methods` must be a list of sem_method() objects")
  benchmarks <- if (is.null(benchmarks)) character(0) else match.arg(benchmarks, c("mean", "lm"), several.ok = TRUE)
  n <- nrow(data)
  if (k < 2 || k > n) stop("`k` must be between 2 and nrow(data)")
  data <- as.data.frame(data)

  # result columns: one per (method, construction), plus benchmarks
  cols <- character(0); col_method <- character(0); col_constr <- character(0)
  for (nm in names(methods)) {
    m <- methods[[nm]]; lab <- if (is.null(m$label)) nm else m$label
    for (cc in m$construction) {
      cols <- c(cols, if (length(m$construction) > 1) paste(lab, cc, sep = "_") else lab)
      col_method <- c(col_method, nm); col_constr <- c(col_constr, cc)
    }
  }
  for (b in benchmarks) { cols <- c(cols, b); col_method <- c(col_method, ""); col_constr <- c(col_constr, b) }
  if (anyDuplicated(cols)) stop("duplicate result column names: ", paste(cols[duplicated(cols)], collapse = ", "))
  nc <- length(cols); ny <- length(ynames)

  if (!is.null(seed)) set.seed(seed)
  folds <- matrix(NA_integer_, n, reps)
  for (r in seq_len(reps)) folds[, r] <- sample(rep(seq_len(k), length.out = n))

  se_item <- array(NA_real_, c(n, reps, nc, ny), dimnames = list(NULL, NULL, cols, ynames))
  yhat <- se_item
  se_obs <- array(NA_real_, c(n, reps, nc), dimnames = list(NULL, NULL, cols))
  mse_in <- array(NA_real_, c(reps, k, nc), dimnames = list(NULL, NULL, cols))
  fail <- stats::setNames(integer(nc), cols)
  reasons <- stats::setNames(vector("list", nc), cols)

  record <- function(col, r, te, out, Yte) {
    out <- as.matrix(out)[, ynames, drop = FALSE]
    se_item[te, r, col, ] <<- (Yte - out)^2
    yhat[te, r, col, ] <<- out
    se_obs[te, r, col] <<- rowMeans((Yte - out)^2)
  }
  note_fail <- function(col, msg) { fail[col] <<- fail[col] + 1L; reasons[[col]] <<- c(reasons[[col]], msg) }

  for (r in seq_len(reps)) {
    for (f in seq_len(k)) {
      if (verbose) cat(sprintf("rep %d fold %d\n", r, f))
      te <- folds[, r] == f
      train <- data[!te, , drop = FALSE]; test <- data[te, , drop = FALSE]
      Ytr <- as.matrix(train[, ynames, drop = FALSE]); Yte <- as.matrix(test[, ynames, drop = FALSE])

      for (nm in names(methods)) {
        m <- methods[[nm]]
        my_cols <- cols[col_method == nm]; my_constr <- col_constr[col_method == nm]
        fit <- tryCatch(m$fit(train), error = function(e) e)
        if (inherits(fit, "error") || is.null(fit)) {
          for (cc in my_cols) note_fail(cc, if (is.null(fit)) "fit returned NULL" else conditionMessage(fit)); next
        }
        if (inherits(fit, "lavaan") && !isTRUE(lavaan::lavInspect(fit, "converged"))) {
          for (cc in my_cols) note_fail(cc, "lavaan did not converge"); next
        }
        par <- tryCatch(do.call(sem_params, c(list(fit), m$params_args)), error = function(e) e)
        if (inherits(par, "error")) { for (cc in my_cols) note_fail(cc, conditionMessage(par)); next }
        for (i in seq_along(my_cols)) {
          col <- my_cols[i]; cc <- my_constr[i]
          out <- tryCatch(predict_oos(par, test, ynames, xnames, construction = cc), error = function(e) e)
          if (inherits(out, "error")) { note_fail(col, conditionMessage(out)); next }
          record(col, r, te, out, Yte)
          ins <- tryCatch(predict_oos(par, train, ynames, xnames, construction = cc), error = function(e) NULL)
          if (!is.null(ins)) mse_in[r, f, col] <- mean((Ytr - as.matrix(ins)[, ynames, drop = FALSE])^2)
        }
      }
      for (b in benchmarks) {
        if (b == "mean") {
          out <- matrix(colMeans(Ytr), nrow(test), ny, byrow = TRUE, dimnames = list(NULL, ynames))
          ins <- matrix(colMeans(Ytr), nrow(train), ny, byrow = TRUE)
        } else {
          fm <- stats::as.formula(sprintf("cbind(%s) ~ %s", paste(ynames, collapse = ","), paste(xnames, collapse = "+")))
          lmfit <- tryCatch(stats::lm(fm, data = train), error = function(e) e)
          if (inherits(lmfit, "error")) { note_fail(b, conditionMessage(lmfit)); next }
          out <- stats::predict(lmfit, newdata = test); ins <- stats::predict(lmfit, newdata = train)
        }
        record(b, r, te, out, Yte)
        mse_in[r, f, b] <- mean((Ytr - as.matrix(ins))^2)
      }
    }
  }
  structure(list(data = data, ynames = ynames, xnames = xnames, k = k, reps = reps, seed = seed,
                 folds = folds, columns = cols, column_method = col_method,
                 column_construction = col_constr, se_item = se_item, se_obs = se_obs,
                 yhat = yhat, mse_in = mse_in, fail = fail, reasons = reasons),
            class = "sem_cv")
}

#' @export
print.sem_cv <- function(x, ...) {
  cat(sprintf("<sem_cv> %d x %d-fold cross-validation, n = %d, seed = %s\n", x$reps, x$k,
              nrow(x$data), if (is.null(x$seed)) "none" else x$seed))
  cat("  predict:", paste(x$ynames, collapse = ", "), "\n  from:   ", paste(x$xnames, collapse = ", "), "\n")
  cat("  columns:", paste(x$columns, collapse = ", "), "\n")
  if (any(x$fail > 0)) {
    cat("  failed folds:", paste(sprintf("%s %d", names(x$fail)[x$fail > 0], x$fail[x$fail > 0]), collapse = "; "), "\n")
  } else cat("  failed folds: none\n")
  invisible(x)
}
