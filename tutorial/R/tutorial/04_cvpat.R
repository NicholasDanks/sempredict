# ==============================================================================
# 04_cvpat.R
# CVPAT (Liengaard et al., 2021) applied as published: ONE k-fold cross-validation,
# per-observation loss = squared error averaged over the outcome indicators,
# one-sided paired t-test of the loss difference benchmark - model over the held-out
# cases (H1: the model's loss is lower than the benchmark's, as in Liengaard et al., 2021).
#
# Reads the saved repeated-CV objects written by 01_ and 02_ (same seed, same folds)
# and reports CVPAT on repetition 1. The other 19 repetitions are reported only as a
# sensitivity range for the p value, because training sets overlap across repetitions
# and the repetitions are therefore not independent evidence.
# Run from the project root:  Rscript R/tutorial/04_cvpat.R
# ==============================================================================
out_dir <- "results/tutorial"
sink(file.path(out_dir, "04_cvpat.txt"), split = TRUE)
cat("CVPAT on a single 10-fold repetition —", format(Sys.time()), "\n\n")

cvpat <- function(se_obs, model, bench, r) {
  d <- se_obs[, r, bench] - se_obs[, r, model]      # benchmark loss - model loss
  d <- d[is.finite(d)]
  tt <- t.test(d, alternative = "greater")   # one-sided: H1 mean(d) > 0
  c(n = length(d), d = mean(d), t = unname(tt$statistic), p = tt$p.value)
}
report <- function(res, label, methods = c("cbsem", "pls_chain", "plsc_chain", "plsc_cov_blp", "lm")) {
  cat(sprintf("--- %s ---\n", label))
  for (bench in c("lm", "mean")) {
    cat(sprintf("Benchmark: %s\n", bench))
    tab <- t(sapply(setdiff(methods, bench), function(m) {
      r1 <- cvpat(res$se_obs, m, bench, 1)
      ps <- sapply(seq_len(dim(res$se_obs)[2]), function(r) cvpat(res$se_obs, m, bench, r)["p"])
      c(round(r1, 3), p_min_20reps = round(min(ps), 3), p_max_20reps = round(max(ps), 3),
        reps_p_lt_05 = sum(ps < .05))
    }))
    print(tab)
  }
  cat("\n")
}
x <- readRDS(file.path(out_dir, "01_political_democracy.rds"))
report(x$results$DA, "PoliticalDemocracy, direct antecedents (y5-y8 from x1-x3, y1-y4)")
report(x$results$EA, "PoliticalDemocracy, earliest antecedent (y5-y8 from x1-x3)")
y <- readRDS(file.path(out_dir, "02_mobi_low_reliability.rds"))
report(y$case1$results[[1]], "mobi case 1: EXP -> SAT")
report(y$case2$results[[1]], "mobi case 2: IMAG + EXP -> SAT")
cat("n = held-out cases with a finite loss for both model and benchmark on repetition 1\n")
cat("(PLSc covariance BLP: cases in inadmissible folds are excluded, so its n is smaller).\n")
cat("d = mean(benchmark loss - model loss); positive favours the model.\n")
cat("p is one-sided (H1: model loss below benchmark loss); p > .5 means d < 0.\n")
sink()
