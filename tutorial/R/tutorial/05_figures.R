# ==============================================================================
# 05_figures.R
# Figure 1 (article): the workflow and the two prediction constructions (schematic).
# Figure 2 (article): observed vs predicted satisfaction items for the two PLSc
#           constructions on the same held-out cases (mobi, EXP -> SAT,
#           repetition 1 of the 20 x 10-fold run saved by 02_mobi_low_reliability.R).
# Run from the project root:  Rscript R/tutorial/05_figures.R
# ==============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(grid) })
dir.create("figures", showWarnings = FALSE)

# ---- Figure 1: schematic -----------------------------------------------------
box <- function(x, y, w, h, label, fill = "white", size = 3.1, col = "grey25", fontface = "plain") {
  list(annotate("rect", xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2,
                fill = fill, colour = col, linewidth = .4),
       annotate("text", x = x, y = y, label = label, size = size, lineheight = .9, fontface = fontface))
}
arrow_ <- function(x1, y1, x2, y2, label = NULL, above = TRUE, size = 2.7) {
  out <- list(annotate("segment", x = x1, y = y1, xend = x2, yend = y2, linewidth = .5,
                       arrow = arrow(length = unit(1.8, "mm"), type = "closed"), colour = "grey25"))
  if (!is.null(label)) out <- c(out, list(annotate("text", x = (x1 + x2)/2, y = (y1 + y2)/2 + if (above) .45 else -.45,
                                                    label = label, size = size, lineheight = .9, colour = "grey15")))
  out
}
p1 <- ggplot() + xlim(0, 20) + ylim(0, 12.5) + theme_void() + coord_fixed(ratio = 1, clip = "off")
# workflow strip
steps <- c("1 Define the\ntarget", "2 Choose the\ndesign", "3 Fit the SEM\n(ML, PLS, PLSc)",
           "4 Turn the fit into\npredictions", "5 Cross-validate\n(same folds)", "6 Compare with\nbenchmarks", "7 Diagnose\nand report")
xs <- seq(1.6, 18.4, length.out = 7)
p1 <- p1 + annotate("text", x = 0.2, y = 12.2, label = "A. The workflow", hjust = 0, fontface = "bold", size = 3.6)
for (i in seq_along(steps)) p1 <- p1 + box(xs[i], 10.9, 2.5, 1.5, steps[i], fill = if (i == 4) "#FFF2CC" else "grey96")
for (i in 1:6) p1 <- p1 + arrow_(xs[i] + 1.28, 10.9, xs[i + 1] - 1.28, 10.9)
# construction 1: score chain
p1 <- p1 + annotate("text", x = 0.2, y = 8.9, label = "B. Step 4, construction 1: the construct-score chain (PLSpredict)", hjust = 0, fontface = "bold", size = 3.6) +
  box(2.4, 7.4, 3.8, 1.7, "New case's\npredictor items", fill = "#DCE9F5") +
  box(7.4, 7.4, 3.2, 1.7, "Predictor\nconstruct score", fill = "grey96") +
  box(12.2, 7.4, 3.2, 1.7, "Predicted outcome\nconstruct score", fill = "grey96") +
  box(17.4, 7.4, 3.8, 1.7, "Predicted\noutcome items", fill = "#DCE9F5") +
  arrow_(4.35, 7.4, 5.75, 7.4, "weights W") + arrow_(9.05, 7.4, 10.55, 7.4, "paths B") + arrow_(13.85, 7.4, 15.45, 7.4, "loadings Λ") +
  annotate("label", x = 9.8, y = 5.6, size = 2.8, lineheight = .95, fill = "#FFF2CC",
           label = "PLSc corrects B and Λ for unreliability but leaves the score uncorrected.\nCorrected paths applied to an uncorrected score over-disperse the prediction:\none predicted point buys fewer observed points (slope < 1, Figure 2, left).")
# construction 2: covariance-based
p1 <- p1 + annotate("text", x = 0.2, y = 3.9, label = "C. Step 4, construction 2: covariance-based prediction (de Rooij et al., 2023)", hjust = 0, fontface = "bold", size = 3.6) +
  box(2.4, 2.4, 3.8, 1.7, "New case's\npredictor items", fill = "#DCE9F5") +
  box(9.8, 2.4, 6.6, 1.7, "Item covariances the fitted model implies\nΣ = ΛΦΛ′ + Θ  (loadings, construct correlations, error)", fill = "grey96", size = 2.9) +
  box(17.4, 2.4, 3.8, 1.7, "Predicted\noutcome items", fill = "#DCE9F5") +
  arrow_(4.35, 2.4, 6.45, 2.4) + arrow_(13.15, 2.4, 15.45, 2.4, "regression implied\nby Σ (Eq. 1)") +
  annotate("label", x = 9.8, y = 0.55, size = 2.8, lineheight = .95, fill = "#FFF2CC",
           label = "For PLSc, Σ is built from the corrected loadings and construct correlations, so scores and parameters match.\nIt exists only when the corrected correlations are admissible (no implied correlation above one); otherwise the fold is reported as unavailable.")
ggsave("figures/fig1_constructions.png", p1, width = 9.2, height = 5.9, dpi = 300, bg = "white")

# ---- Figure 2: observed vs predicted, mobi case 1, repetition 1 --------------
y <- readRDS("results/tutorial/02_mobi_low_reliability.rds"); r <- y$case1$results[[1]]
d <- seminr::mobi; Y <- as.matrix(d[, r$y])
mk <- function(m, lab) {
  yh <- r$yhat[, 1, m, ]
  yc <- scale(Y, scale = FALSE); hc <- scale(yh, scale = FALSE)
  data.frame(construction = lab, observed = as.vector(yc), predicted = as.vector(hc),
             slope = sum(yc * hc) / sum(hc^2))
}
dd <- rbind(mk("plsc_chain", "PLSc, construct-score chain"), mk("plsc_cov_blp", "PLSc, covariance-based"))
dd$construction <- factor(dd$construction, levels = unique(dd$construction))
lab <- unique(dd[, c("construction", "slope")])
lab$text <- sprintf("slope = %.2f\none predicted point ≈ %.2f observed points", lab$slope, lab$slope)
p2 <- ggplot(dd, aes(predicted, observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(alpha = .35, size = 1.2, colour = "grey30") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#C0392B", linewidth = .9) +
  geom_text(data = lab, aes(x = -2.9, y = 3.7, label = text), hjust = 0, size = 3, lineheight = .95, colour = "#C0392B") +
  facet_wrap(~ construction) + coord_equal(xlim = c(-3.2, 3.2), ylim = c(-4.2, 4.2)) +
  labs(x = "Predicted satisfaction item (difference from the item mean, points)",
       y = "Observed satisfaction item\n(difference from the item mean, points)") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
ggsave("figures/fig2_dispersion.png", p2, width = 8, height = 4.6, dpi = 300, bg = "white")
cat("slopes drawn:", paste(sprintf("%s %.3f", lab$construction, lab$slope), collapse = "; "), "\n")
# Journal-ready TIFF copies (300 dpi) via ImageMagick when available.
for (f in c("fig1_constructions", "fig2_dispersion"))
  if (nzchar(Sys.which("magick"))) system2("magick", c(sprintf("figures/%s.png", f), "-units", "PixelsPerInch", "-density", "300", sprintf("figures/%s.tiff", f)))
cat("Written figures/fig1_constructions.{png,tiff}, figures/fig2_dispersion.{png,tiff}\n")
