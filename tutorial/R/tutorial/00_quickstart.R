# ==============================================================================
# 00_quickstart.R
# Minimal reader-facing example for the Teacher's Corner tutorial.
#
# This script shows how to fit the PoliticalDemocracy model with
#   1) covariance-based SEM using lavaan,
#   2) PLS using seminr, and
#   3) consistent PLS using seminr,
# and how to generate the prediction constructions used in the article.
#
# Run from the project root after installing lavaan and seminr.
# The full repeated 20 x 10-fold analysis is in 01_political_democracy.R.
# ==============================================================================

suppressPackageStartupMessages({
  library(lavaan)
  library(seminr)
})

source("R/tutorial/00_helpers.R")
set.seed(20260825)

data <- lavaan::PoliticalDemocracy

blocks <- list(
  ind60 = paste0("x", 1:3),
  dem60 = paste0("y", 1:4),
  dem65 = paste0("y", 5:8)
)
causal_order <- c("ind60", "dem60", "dem65")
exo_constructs <- "ind60"
y_names <- blocks$dem65
x_DA <- c(blocks$ind60, blocks$dem60)  # direct-antecedent design

lav_syntax <- '
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
  dem60 ~ ind60
  dem65 ~ ind60 + dem60
'

sm <- relationships(
  paths(from = "ind60", to = "dem60"),
  paths(from = c("ind60", "dem60"), to = "dem65")
)

mm_pls <- constructs(
  composite("ind60", blocks$ind60, weights = mode_A),
  composite("dem60", blocks$dem60, weights = mode_A),
  composite("dem65", blocks$dem65, weights = mode_A)
)

mm_plsc <- constructs(
  reflective("ind60", blocks$ind60),
  reflective("dem60", blocks$dem60),
  reflective("dem65", blocks$dem65)
)

# A simple train/test split for illustration only.
# The manuscript results use repeated 10-fold cross-validation.
id <- sample(seq_len(nrow(data)), 60)
train <- data[id, ]
test  <- data[-id, ]

# Fit the three SEM estimators on the training data.
fit_ml   <- sem(lav_syntax, data = train, meanstructure = TRUE)
fit_pls  <- estimate_pls(train, mm_pls,  sm)
fit_plsc <- estimate_pls(train, mm_plsc, sm)

# 1) Covariance-based best linear prediction from the fitted CB-SEM.
yhat_ml <- lavPredictY(
  fit_ml,
  newdata = test,
  ynames = y_names,
  xnames = x_DA
)

# 2) Native PLS construct-score chain.
p_pls <- seminr_params(fit_pls, FALSE, blocks, causal_order, exo_constructs)
yhat_pls <- predict_chain(
  p_pls$W, p_pls$B, p_pls$L,
  test, train, x_DA, y_names, blocks, causal_order
)

# 3) PLSc construct-score chain and, when admissible,
#    covariance-based prediction from the fitted common-factor solution.
p_plsc <- seminr_params(fit_plsc, TRUE, blocks, causal_order, exo_constructs)
yhat_plsc_chain <- predict_chain(
  p_plsc$W, p_plsc$B, p_plsc$L,
  test, train, x_DA, y_names, blocks, causal_order
)

if (params_ok(p_plsc)) {
  yhat_plsc_cov <- predict_cov_blp(
    p_plsc$L, p_plsc$Phi, p_plsc$Theta,
    test, train, x_DA, y_names
  )
} else {
  yhat_plsc_cov <- NULL
  message("PLSc covariance-based prediction unavailable: ", p_plsc$inadmissible)
}

# Simple RMSE summaries for the held-out cases.
rmse <- function(y, yhat) sqrt(mean((as.matrix(y) - as.matrix(yhat))^2))
Ytest <- test[, y_names]

cat(sprintf("ML covariance BLP RMSE: %.3f\n", rmse(Ytest, yhat_ml)))
cat(sprintf("PLS score-chain RMSE:    %.3f\n", rmse(Ytest, yhat_pls)))
cat(sprintf("PLSc score-chain RMSE:   %.3f\n", rmse(Ytest, yhat_plsc_chain)))
if (!is.null(yhat_plsc_cov)) {
  cat(sprintf("PLSc covariance BLP RMSE: %.3f\n", rmse(Ytest, yhat_plsc_cov)))
}

cat("\nFor the repeated cross-validation, benchmarks, OFR, dispersion diagnostic,\n")
cat("and paired loss comparison used in the article, run:\n")
cat("  Rscript R/tutorial/01_political_democracy.R\n")
