# Evaluation: discrimination (AUC), pairwise comparisons and ordinal trend

#' Evaluate scores against sample labels
#'
#' For each score: (i) AUC of NT versus all non-NT samples, (ii) pairwise comparisons between adjacent
#' and all ordered groups (Wilcoxon rank-sum, AUC with DeLong 95% CI, direction fixed so that a higher
#' stage has a higher score; Benjamini-Hochberg within score), and (iii) for three or more ordered
#' groups a Spearman correlation and a linear-by-linear trend test (`coin`) with the stage.
#' Groups with <= `min_n` samples are skipped (reported as NA).
#'
#' @param samples data.frame, e.g. `hixscore()$samples`, with a group column and score columns.
#' @param scores Score columns to evaluate.
#' @param group_col Column holding the group label.
#' @param stages Ordered character vector of groups from least to most advanced, e.g.
#'   `c("NT", "non-malignant", "malignant")`; or a named vector mapping raw labels to those
#'   (e.g. `c(NT = "NT", OLP = "non-malignant", Tumor = "malignant")`).
#' @param min_n Groups with this many samples or fewer are skipped.
#' @param exclude_flagged If `TRUE` and `samples` has the column `M_epi_below_training_min`, samples flagged
#'   there (epithelial cell number below the training range, see [hixscore()]) are left out of the analysis
#'   of **all** scores.
#' @return list(auc_nt_vs_disease, pairwise, trend).
#' @export
evaluate_scores <- function(samples, scores = c("E_seHIX_z", "S_C_z", "S_M_z", "sysHIX"),
                            group_col = "tissue", stages = c("NT", "non-malignant", "malignant"),
                            min_n = 3, exclude_flagged = FALSE) {
  .need("pROC", "evaluate_scores()")
  if (exclude_flagged && "M_epi_below_training_min" %in% names(samples))
    samples <- samples[!(samples$M_epi_below_training_min %in% TRUE), , drop = FALSE]
  raw <- as.character(samples[[group_col]])
  stage <- if (!is.null(names(stages))) unname(stages[raw]) else raw
  lv <- if (!is.null(names(stages))) unique(unname(stages)) else stages
  d <- samples[!is.na(stage) & stage %in% lv, , drop = FALSE]; d$stage <- factor(stage[!is.na(stage) & stage %in% lv], lv, ordered = TRUE)
  scores <- intersect(scores, names(d))
  nstage <- as.integer(d$stage)
  one_auc <- function(x, hi, lo) {   # AUC that group hi scores higher than lo
    k <- !is.na(x)
    a <- d$stage[k] %in% c(lo, hi); xx <- x[k][a]; yy <- as.integer(d$stage[k][a] %in% hi)
    if (sum(yy == 1) <= min_n || sum(yy == 0) <= min_n) return(c(n_lo = sum(yy == 0), n_hi = sum(yy == 1), auc = NA, lo = NA, hi = NA, p = NA))
    r <- pROC::roc(yy, xx, direction = "<", quiet = TRUE); ci <- suppressWarnings(as.numeric(pROC::ci.auc(r, method = "delong")))
    c(n_lo = sum(yy == 0), n_hi = sum(yy == 1), auc = as.numeric(pROC::auc(r)), lo = ci[1], hi = ci[3],
      p = stats::wilcox.test(xx[yy == 1], xx[yy == 0], exact = FALSE)$p.value)
  }
  auc_tab <- do.call(rbind, lapply(scores, function(s) data.frame(score = s, t(one_auc(d[[s]], lv[-1], lv[1])))))
  pairs <- utils::combn(seq_along(lv), 2, simplify = FALSE)
  pw <- do.call(rbind, lapply(scores, function(s) {
    t <- do.call(rbind, lapply(pairs, function(p) data.frame(score = s, lower = lv[p[1]], higher = lv[p[2]], t(one_auc(d[[s]], lv[p[2]], lv[p[1]])))))
    t$p_adj_BH <- stats::p.adjust(t$p, "BH"); t }))
  tr <- NULL
  if (length(lv) >= 3) tr <- do.call(rbind, lapply(scores, function(s) {
    k <- !is.na(d[[s]]); x <- d[[s]][k]; g <- d$stage[k]
    ok <- all(table(g) > min_n)
    p <- if (ok && requireNamespace("coin", quietly = TRUE))
      as.numeric(coin::pvalue(coin::independence_test(x ~ as.numeric(g), distribution = "asymptotic"))) else NA
    data.frame(score = s, n = sum(k), spearman_rho = if (ok) unname(stats::cor(x, as.numeric(g), method = "spearman")) else NA,
               trend_p = p)
  }))
  list(auc_nt_vs_disease = auc_tab, pairwise = pw, trend = tr)
}

#' Boxplot of scores by group with sample numbers
#' @param samples,scores,group_col,stages as in [evaluate_scores()].
#' @return A ggplot object.
#' @export
plot_scores <- function(samples, scores = c("E_seHIX_z", "S_C_z", "S_M_z", "sysHIX"),
                        group_col = "tissue", stages = NULL) {
  .need("ggplot2", "plot_scores()")
  scores <- intersect(scores, names(samples))
  g <- as.character(samples[[group_col]]); if (!is.null(stages)) g <- ifelse(g %in% names(stages), unname(stages[g]), NA)
  lv <- if (is.null(stages)) unique(g[!is.na(g)]) else unique(unname(stages))
  long <- do.call(rbind, lapply(scores, function(s) data.frame(group = g, score = s, value = samples[[s]])))
  long <- long[!is.na(long$group) & !is.na(long$value), ]
  n <- table(long$score, long$group)
  long$label <- factor(paste0(long$group, "\n(n=", n[cbind(long$score, long$group)], ")"))
  long$label <- factor(long$label, levels = unique(long$label[order(match(long$group, lv))]))
  ggplot2::ggplot(long, ggplot2::aes(label, value, fill = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
    ggplot2::facet_wrap(~score, scales = "free", nrow = 1) + ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = "Score", fill = NULL) + ggplot2::theme(legend.position = "none")
}
