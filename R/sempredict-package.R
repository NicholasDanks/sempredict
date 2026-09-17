#' sempredict: out-of-sample prediction assessment for SEM
#'
#' One assessment layer over several estimators. Fit the model with lavaan
#' (ML), seminr (PLS) or seminr with `reflective()` constructs (consistent
#' PLS); [sem_params()] extracts what prediction needs; [predict_oos()] turns
#' the fit into case-level predictions by the construct-score chain or by the
#' model-implied regression; [cv_predict()] runs repeated k-fold
#' cross-validation with identical folds for every method; [assess()],
#' [cvpat()], [failures()] and [plot_dispersion()] report the result.
#'
#' @keywords internal
#' @importFrom stats sd cor setNames complete.cases
#' @importFrom utils globalVariables
"_PACKAGE"

utils::globalVariables(c("predicted", "observed", "text"))
