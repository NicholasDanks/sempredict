# ==============================================================================
# 02_mobi_low_reliability.R
# Metric alignment in PLSc score-based prediction at a magnitude that matters: seminr's `mobi`
# (ECSI mobile-phone customer satisfaction, n = 250).
#
# Case 1  single antecedent:        EXP (CUEX1-3) -> SAT (CUSA1-3)
#         rho_A(EXP) ~ .46, so the PLSc chain should be over-dispersed by
#         1/sqrt(rho_A) ~ 1.47 relative to the best linear prediction based on the same composite score.
# Case 2  two correlated antecedents of unequal reliability:
#         IMAG (IMAG1-5) + EXP (CUEX1-3) -> SAT
#         here the mismatch is a REWEIGHTING: the chain changes the relative weighting of antecedent blocks when their reliabilities differ.
#
# Same machinery as 01_ (source 00_helpers.R).  Run from the project root.
# ==============================================================================

suppressPackageStartupMessages({ library(lavaan); library(seminr) })
source("R/tutorial/00_helpers.R")

SEED <- 20260825; K_FOLDS <- 10; N_REPS <- 20
out_dir <- "results/tutorial"; if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
sink(file.path(out_dir, "02_mobi_low_reliability.txt"), split = TRUE)
cat("mobi low-reliability run —", format(Sys.time()), "\n")
cat("lavaan", as.character(packageVersion("lavaan")), "| seminr", as.character(packageVersion("seminr")), "\n\n")

data <- seminr::mobi; n <- nrow(data)
methods <- c("cbsem", "pls_chain", "pls_cov_plugin", "plsc_chain", "plsc_cov_blp", "lm", "mean")

run_case <- function(label, blocks, exo, order, lav_syntax, sm, designs) {
  cat(sprintf("\n==================== %s ====================\n", label))
  assign("blocks", blocks, envir = .GlobalEnv); assign("exo_constructs", exo, envir = .GlobalEnv)
  assign("causal_order", order, envir = .GlobalEnv); assign("lav_syntax", lav_syntax, envir = .GlobalEnv)
  assign("sm", sm, envir = .GlobalEnv); assign("designs", designs, envir = .GlobalEnv)
  assign("mm_pls",  do.call(constructs, lapply(names(blocks), function(cn) composite(cn, blocks[[cn]], weights = mode_A))), envir = .GlobalEnv)
  assign("mm_plsc", do.call(constructs, lapply(names(blocks), function(cn) reflective(cn, blocks[[cn]]))), envir = .GlobalEnv)

  # --- full-sample mechanism table
  m_pls <- estimate_pls(data, mm_pls, sm); m_plsc <- estimate_pls(data, mm_plsc, sm)
  rA <- sapply(order, function(cn) seminr::rho_A(m_plsc, cn))
  cat("rho_A:", paste(sprintf("%s = %.3f", order, rA), collapse = "; "), "\n")
  cat("1/sqrt(rho_A):", paste(sprintf("%s = %.3f", order, 1 / sqrt(rA)), collapse = "; "), "\n")
  x_names <- designs[[1]]$x; y_names <- designs[[1]]$y
  pp <- seminr_params(m_pls, FALSE, blocks, order, exo)
  pc <- tryCatch(seminr_params(m_plsc, TRUE, blocks, order, exo), error = function(e) e)
  fit <- sem(lav_syntax, data = data, meanstructure = TRUE)
  imp <- fitted(fit); G_ml <- imp$cov[y_names, x_names] %*% solve(imp$cov[x_names, x_names])
  # put ML Gamma on the standardised-x / standardised-y scale for comparability
  sx <- apply(data[, x_names], 2, sd); sy <- apply(data[, y_names], 2, sd)
  G_ml_z <- sweep(sweep(G_ml, 2, sx, "*"), 1, sy, "/")
  cat("\nBlockwise prediction-coefficient magnitude (mean over SAT items of sum|Gamma| within block):\n")
  tab <- rbind(ML_rule     = block_credit(G_ml_z, blocks, x_names),
               PLS_chain   = block_credit(gamma_chain(pp, x_names, y_names, blocks, order), blocks, x_names),
               PLS_cov_plugin    = block_credit(gamma_sem(pp, x_names, y_names), blocks, x_names))
  if (params_ok(pc)) tab <- rbind(tab,
               PLSc_chain  = block_credit(gamma_chain(pc, x_names, y_names, blocks, order), blocks, x_names),
               PLSc_cov_BLP   = block_credit(gamma_sem(pc, x_names, y_names), blocks, x_names)) else
    cat("PLSc full-sample fit inadmissible:", if (inherits(pc, "error")) conditionMessage(pc) else pc$inadmissible, "\n")
  print(round(tab, 3))
  if (params_ok(pc)) {
    cat("PLSc_chain / PLSc_cov_BLP coefficient-magnitude ratio by block:",
        paste(sprintf("%s %.3f", colnames(tab), tab["PLSc_chain", ] / tab["PLSc_cov_BLP", ]), collapse = "; "), "\n")
    cat("PLSc_chain / ML_rule coefficient-magnitude ratio by block:  ",
        paste(sprintf("%s %.3f", colnames(tab), tab["PLSc_chain", ] / tab["ML_rule", ]), collapse = "; "), "\n")
    cat("SD of full-sample predictions, PLSc_chain / PLSc_cov_BLP (per SAT item):",
        paste(round(apply(predict_chain(pc$W, pc$B, pc$L, data, data, x_names, y_names, blocks, order), 2, sd) /
                    apply(predict_cov_blp(pc$L, pc$Phi, pc$Theta, data, data, x_names, y_names), 2, sd), 3), collapse = " "), "\n")
  }
  # --- CV
  cat(sprintf("\n%d x %d-fold CV, seed %d\n", N_REPS, K_FOLDS, SEED))
  res <- lapply(names(designs), run_design); names(res) <- names(designs)
  tabs <- lapply(res, summarise_design)
  list(results = res, tables = tabs, credit = tab, rhoA = rA)
}

# ---- Case 1 ------------------------------------------------------------------
case1 <- run_case(
  "Case 1: EXP -> SAT (single antecedent, low reliability)",
  blocks = list(EXP = paste0("CUEX", 1:3), SAT = paste0("CUSA", 1:3)),
  exo = "EXP", order = c("EXP", "SAT"),
  lav_syntax = 'EXP =~ CUEX1 + CUEX2 + CUEX3\n SAT =~ CUSA1 + CUSA2 + CUSA3\n SAT ~ EXP',
  sm = relationships(paths(from = "EXP", to = "SAT")),
  designs = list(DA = list(x = paste0("CUEX", 1:3), y = paste0("CUSA", 1:3))))

# ---- Case 2 ------------------------------------------------------------------
case2 <- run_case(
  "Case 2: IMAG + EXP -> SAT (correlated antecedents, unequal reliability)",
  blocks = list(IMAG = paste0("IMAG", 1:5), EXP = paste0("CUEX", 1:3), SAT = paste0("CUSA", 1:3)),
  exo = c("IMAG", "EXP"), order = c("IMAG", "EXP", "SAT"),
  lav_syntax = 'IMAG =~ IMAG1 + IMAG2 + IMAG3 + IMAG4 + IMAG5\n EXP =~ CUEX1 + CUEX2 + CUEX3\n SAT =~ CUSA1 + CUSA2 + CUSA3\n SAT ~ IMAG + EXP',
  sm = relationships(paths(from = c("IMAG", "EXP"), to = "SAT")),
  designs = list(DA = list(x = c(paste0("IMAG", 1:5), paste0("CUEX", 1:3)), y = paste0("CUSA", 1:3))))

sink()
saveRDS(list(case1 = case1, case2 = case2, seed = SEED), file.path(out_dir, "02_mobi_low_reliability.rds"))
cat("\nWritten to", file.path(out_dir, "02_mobi_low_reliability.txt"), "\n")
