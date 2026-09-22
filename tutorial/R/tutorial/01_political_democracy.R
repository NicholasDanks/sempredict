# ==============================================================================
# 01_political_democracy.R
# Tutorial worked case: out-of-sample prediction assessment for an SEM,
# on lavaan's PoliticalDemocracy (Bollen 1989), n = 75.
#
#   ind60 =~ x1 + x2 + x3
#   dem60 =~ y1 + y2 + y3 + y4
#   dem65 =~ y5 + y6 + y7 + y8
#   dem60 ~ ind60
#   dem65 ~ ind60 + dem60
#
# Two prediction designs (PLSpredict vocabulary):
#   DA  direct antecedents:   predict y5-y8 from x1-x3 AND y1-y4
#   EA  earliest antecedent:  predict y5-y8 from x1-x3 only (through the chain)
#
# Prediction constructions used in the tutorial:
#   chain            composite-score chain W -> B -> Lambda (PLSpredict)
#   covariance BLP   best linear prediction from fitted moments (de Rooij et al. 2023)
#
# Estimators: CB-SEM (ML, covariance BLP via lavaan::lavPredictY),
#             PLS Mode A (native chain; exploratory covariance plug-in retained for diagnostics),
#             PLSc (chain and covariance BLP when admissible), LM, mean.
#
# Assessment: repeated k-fold CV; item RMSE; descriptive paired loss differences; overfit ratio OFR = (MSE_out - MSE_in) / MSE_in; admissibility
# counts for PLSc.  Everything regenerates from this script with the seed below.
#
# Run from the project root:  Rscript R/tutorial/01_political_democracy.R
# ==============================================================================

suppressPackageStartupMessages({
  library(lavaan)
  library(seminr)
})

SEED    <- 20260825
K_FOLDS <- 10
N_REPS  <- 20

out_dir <- "results/tutorial"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

data <- lavaan::PoliticalDemocracy
n <- nrow(data)

# ---- Model -------------------------------------------------------------------
blocks <- list(ind60 = paste0("x", 1:3),
               dem60 = paste0("y", 1:4),
               dem65 = paste0("y", 5:8))
exo_constructs <- "ind60"
causal_order   <- c("ind60", "dem60", "dem65")   # recursive order

lav_syntax <- '
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
  dem60 ~ ind60
  dem65 ~ ind60 + dem60
'
# Note: Bollen's canonical model adds residual covariances (y1~~y5 etc.). They are
# omitted so that the same model is estimable by every method; the CB-SEM fit
# with them is reported once in Section A below for reference.

sm <- relationships(paths(from = "ind60", to = "dem60"),
                    paths(from = c("ind60", "dem60"), to = "dem65"))
mm_pls  <- constructs(composite("ind60", blocks$ind60, weights = mode_A),
                      composite("dem60", blocks$dem60, weights = mode_A),
                      composite("dem65", blocks$dem65, weights = mode_A))
mm_plsc <- constructs(reflective("ind60", blocks$ind60),
                      reflective("dem60", blocks$dem60),
                      reflective("dem65", blocks$dem65))

designs <- list(
  DA = list(x = c(blocks$ind60, blocks$dem60), y = blocks$dem65),
  EA = list(x = blocks$ind60,                  y = blocks$dem65)
)

source("R/tutorial/00_helpers.R")

sink(file.path(out_dir, "01_political_democracy.txt"), split = TRUE)
cat("PoliticalDemocracy tutorial run —", format(Sys.time()), "\n")
cat("lavaan", as.character(packageVersion("lavaan")), "| seminr", as.character(packageVersion("seminr")), "\n\n")

# ---- Section A: one full-sample fit of everything (for the tutorial text) ---
cat("=== A. Full-sample fits ===\n")
fit_full <- sem(lav_syntax, data = data, meanstructure = TRUE)
fit_bollen <- sem(paste(lav_syntax, 'y1 ~~ y5\n y2 ~~ y4 + y6\n y3 ~~ y7\n y4 ~~ y8\n y6 ~~ y8'),
                  data = data, meanstructure = TRUE)
cat(sprintf("CB-SEM (no residual cov): chi2(%d) = %.2f, p = %.3f, CFI = %.3f, RMSEA = %.3f\n",
            fitMeasures(fit_full, "df"), fitMeasures(fit_full, "chisq"),
            fitMeasures(fit_full, "pvalue"), fitMeasures(fit_full, "cfi"),
            fitMeasures(fit_full, "rmsea")))
cat(sprintf("CB-SEM (Bollen, resid cov): chi2(%d) = %.2f, p = %.3f, CFI = %.3f\n",
            fitMeasures(fit_bollen, "df"), fitMeasures(fit_bollen, "chisq"),
            fitMeasures(fit_bollen, "pvalue"), fitMeasures(fit_bollen, "cfi")))
cat("ML latent correlations (boundary proximity: dem60-dem65):\n"); print(round(lavInspect(fit_full, "cor.lv"), 3))
pe <- parameterEstimates(fit_full, standardized = TRUE)
cat("ML standardised paths:\n"); print(pe[pe$op == "~", c("lhs", "rhs", "std.all")], row.names = FALSE)

m_pls  <- estimate_pls(data, mm_pls,  sm)
m_plsc <- estimate_pls(data, mm_plsc, sm)
cat("PLS  paths:\n"); print(round(m_pls$path_coef[c("ind60", "dem60"), c("dem60", "dem65")], 3))
cat("PLSc paths:\n"); print(round(m_plsc$path_coef[c("ind60", "dem60"), c("dem60", "dem65")], 3))
rA_full <- sapply(causal_order, function(cn) seminr::rho_A(m_plsc, cn))
cat("rho_A (full sample):", paste(sprintf("%s = %.3f", causal_order, rA_full), collapse = "; "), "\n")
cat("Predicted 1/sqrt(rho_A) inflation of the PLSc chain, single-antecedent case:",
    sprintf("ind60 %.3f, dem60 %.3f", 1 / sqrt(rA_full["ind60"]), 1 / sqrt(rA_full["dem60"])), "\n")

# Cross-check: lavPredictY (lavaan) == our own implementation of the rule
xd <- designs$DA
own <- {
  imp <- fitted(fit_full)
  Sxx <- imp$cov[xd$x, xd$x]; Syx <- imp$cov[xd$y, xd$x]
  G <- Syx %*% solve(Sxx)
  sweep(as.matrix(data[, xd$x]) %*% t(G), 2, as.numeric(imp$mean[xd$y] - G %*% imp$mean[xd$x]), "+")
}
lav <- lavPredictY(fit_full, newdata = data, ynames = xd$y, xnames = xd$x)
cat(sprintf("lavPredictY vs own de Rooij implementation: max |diff| = %.2e\n\n",
            max(abs(own - lav))))

# ---- Blockwise prediction-coefficient magnitude (DA) -----------------------------
# Under DA, dem65 has TWO correlated antecedents of unequal reliability, so by
# the two-antecedent argument (supplement S2.3) the chain's distortion is a re-weighting rather than a common
# scale factor. Sum |Gamma| within each antecedent block, averaged over the
# outcome items, makes the re-weighting directly visible.
xd_da <- designs$DA
pp <- seminr_params(m_pls,  FALSE, blocks, causal_order, exo_constructs)
pc <- seminr_params(m_plsc, TRUE,  blocks, causal_order, exo_constructs)
imp_da  <- fitted(fit_full)
G_ml    <- imp_da$cov[xd_da$y, xd_da$x] %*% solve(imp_da$cov[xd_da$x, xd_da$x])
sx <- apply(data[, xd_da$x], 2, sd); sy <- apply(data[, xd_da$y], 2, sd)
G_ml_z  <- sweep(sweep(G_ml, 2, sx, "*"), 1, sy, "/")
cred <- rbind(
  ML_rule    = block_credit(G_ml_z, blocks, xd_da$x),
  PLS_chain  = block_credit(gamma_chain(pp, xd_da$x, xd_da$y, blocks, causal_order), blocks, xd_da$x),
  PLS_cov_plugin   = block_credit(gamma_sem(pp, xd_da$x, xd_da$y), blocks, xd_da$x))
if (params_ok(pc)) cred <- rbind(cred,
  PLSc_chain = block_credit(gamma_chain(pc, xd_da$x, xd_da$y, blocks, causal_order), blocks, xd_da$x),
  PLSc_cov_BLP  = block_credit(gamma_sem(pc, xd_da$x, xd_da$y), blocks, xd_da$x))
cat("\nDA blockwise prediction-coefficient magnitude (mean over dem65 items of sum|Gamma|):\n")
print(round(cred, 3))
if (params_ok(pc)) {
  cat("PLSc_chain / PLSc_cov_BLP coefficient-magnitude ratio by block:",
      paste(sprintf("%s %.3f", colnames(cred), cred["PLSc_chain", ] / cred["PLSc_cov_BLP", ]), collapse = "; "), "\n")
  cat("PLSc_chain / ML_rule   coefficient-magnitude ratio by block:",
      paste(sprintf("%s %.3f", colnames(cred), cred["PLSc_chain", ] / cred["ML_rule", ]), collapse = "; "), "\n")
  cat(sprintf("Composite score correlation cor(ind60, dem60) = %.3f; PLSc-disattenuated = %.3f\n",
              cor(m_plsc$construct_scores[, "ind60"], m_plsc$construct_scores[, "dem60"]),
              pc$Phi["ind60", "dem60"]))
  cat("SD of full-sample DA predictions, PLSc_chain / PLSc_cov_BLP (per dem65 item):",
      paste(round(apply(predict_chain(pc$W, pc$B, pc$L, data, data, xd_da$x, xd_da$y, blocks, causal_order), 2, sd) /
                  apply(predict_cov_blp(pc$L, pc$Phi, pc$Theta, data, data, xd_da$x, xd_da$y), 2, sd), 3),
            collapse = " "), "\n")
}

cat(sprintf("=== B. %d x %d-fold CV, seed %d ===\n", N_REPS, K_FOLDS, SEED))
res <- lapply(names(designs), run_design); names(res) <- names(designs)
tabs <- lapply(res, summarise_design)
cat("\n=== C. Common-fold summaries used for Table 2 ===\n")
common_tabs <- lapply(res, summarise_common_plsc_folds)
print(common_tabs)

# The covariance-based prediction is undefined on the inadmissible folds while the score chain
# remains defined. Table 2 therefore uses common_tabs, which rescores every reported method
# and benchmark on the folds where the PLSc covariance BLP is available.

# ---- D. What omitting Bollen's residual covariances costs CB-SEM ---------------------
# The article's specification drops Bollen's (1989) residual covariances so that all three
# estimators fit the same model. Four of them (y1~~y5, y2~~y6, y3~~y7, y4~~y8) link a
# predictor item to an outcome item under the direct-antecedent design, i.e. they sit in
# Sigma_yx. Same folds as above (same seed and draw order), all 200 folds, ML only.
cat("\n=== D. ML with and without Bollen's residual covariances, same folds ===\n")
y_names <- designs$DA$y
lav_bollen <- paste(lav_syntax, "y1 ~~ y5\n y2 ~~ y4 + y6\n y3 ~~ y7\n y4 ~~ y8\n y6 ~~ y8", sep = "\n")
se_b <- array(NA_real_, c(n, N_REPS, 2, 2, length(y_names)),
              dimnames = list(NULL, NULL, c("article", "bollen"), names(designs), y_names))
set.seed(SEED)
for (r in seq_len(N_REPS)) { folds <- sample(rep(seq_len(K_FOLDS), length.out = n))
  for (k in seq_len(K_FOLDS)) { te <- folds == k
    for (m in c("article", "bollen")) {
      f <- sem(if (m == "article") lav_syntax else lav_bollen, data = data[!te, ], meanstructure = TRUE, warn = FALSE)
      if (!lavInspect(f, "converged")) next
      for (d in names(designs)) { p <- lavPredictY(f, newdata = data[te, ], ynames = y_names, xnames = designs[[d]]$x)
        se_b[te, r, m, d, ] <- (as.matrix(data[te, y_names]) - p[, y_names])^2 } } } }
for (d in names(designs)) for (m in c("article", "bollen")) {
  rr <- sapply(seq_len(N_REPS), function(r) mean(sapply(seq_along(y_names), function(j) sqrt(mean(se_b[, r, m, d, j], na.rm = TRUE)))))
  cat(sprintf("%s  %-8s RMSE %.4f (MCSE %.4f); cells without a converged fit: %d\n", d, m, mean(rr), sd(rr) / sqrt(N_REPS),
              sum(is.na(se_b[, , m, d, 1])))) }
sink()
saveRDS(list(results = res, tables = tabs, common_tables = common_tabs, seed = SEED, k = K_FOLDS, reps = N_REPS),
        file.path(out_dir, "01_political_democracy.rds"))
cat("\nWritten to", file.path(out_dir, "01_political_democracy.txt"), "\n")
