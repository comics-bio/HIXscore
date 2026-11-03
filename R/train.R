# Mode 3: train the layer models and the fusion model on user data

#' Assemble the sample-level training table
#'
#' Joins E, C and M layer outputs with sample labels. Samples need all three layers.
#' @param E,C,M Outputs of [compute_E_seHIX()], [compute_S_C()] (needs the C feature columns)
#'   and [compute_M_CM_scores()].
#' @param labels data.frame with `sample_id`, `label` and (recommended) `series`.
#' @export
build_training_table <- function(E, C, M, labels) {
  stopifnot(all(c("sample_id", "label") %in% names(labels)))
  cm <- grep("^M_CM\\d+$", names(M), value = TRUE)
  d <- merge(C[, c("sample_id", "n_epi_cells", "C_Low_Frac", "C_Mid_Frac", "C_High_Frac", "C_mean")],
             E[, c("sample_id", "E_seHIX")], by = "sample_id")
  d <- merge(d, M[, c("sample_id", cm)], by = "sample_id")
  merge(d, labels[, intersect(c("sample_id", "label", "series"), names(labels))], by = "sample_id")
}

.grouped_folds <- function(y, groups, nFolds, seed) {
  set.seed(seed)
  if (!is.null(groups) && length(unique(groups)) >= nFolds) {
    g <- unique(groups); g <- g[sample(length(g))]
    fold_of <- stats::setNames(rep_len(seq_len(nFolds), length(g)), g)
    f <- fold_of[as.character(groups)]
    lapply(seq_len(nFolds), function(k) which(f == k))
  } else {
    if (!is.null(groups)) warning("Fewer series than folds; falling back to stratified folds (series leakage possible).", call. = FALSE)
    unname(caret::createFolds(y, k = nFolds, list = TRUE))
  }
}

.fit_layer <- function(X, y, nFolds, alpha) {
  tune <- ordinalNet::ordinalNetTune(X, y, family = "cumulative", link = "logit", alpha = alpha,
                                     nFolds = nFolds, standardize = FALSE, printProgress = FALSE)
  lam <- tune$lambdaVals[which.max(rowMeans(tune$loglik))]
  b <- stats::coef(tune$fit, matrix = TRUE, lambda = lam)
  list(fit = tune$fit, lambda = lam, beta = b[rownames(b) != "(Intercept)", 1])
}

.oof_score <- function(X, y, lambda, folds, alpha) {
  oof <- rep(NA_real_, nrow(X))
  for (h in folds) {
    tr <- setdiff(seq_len(nrow(X)), h)
    f <- ordinalNet::ordinalNet(X[tr, , drop = FALSE], y[tr], family = "cumulative", link = "logit",
                                alpha = alpha, lambdaVals = lambda, standardize = FALSE)
    b <- stats::coef(f, matrix = TRUE); b <- b[rownames(b) != "(Intercept)", 1]
    oof[h] <- -as.numeric(X[h, , drop = FALSE] %*% b)
  }
  oof
}

#' Train a sysHIX model
#'
#' Two-stage ordinal-logistic model:
#' (1) regress out log(epithelial cell number) with coefficients fitted on the training set,
#' (2) fit the C-layer and M-layer models on z-scored features (z-scored with NT samples),
#' (3) get **out-of-fold** S_C / S_M for the training samples (folds grouped by `series`),
#' (4) fit the fusion model on E, S_C, S_M. sysHIX is the expected class index.
#'
#' With three ordered classes (e.g. NT < non-malignant < malignant) sysHIX lies in 0 to 2;
#' with two classes (NT vs disease) the fusion is a penalised logistic regression and sysHIX
#' is the probability of disease in 0 to 1.
#'
#' @param tab Table from [build_training_table()].
#' @param levels Ordered class labels in `tab$label`; length 2 or 3, first = NT.
#' @param scHIX_thresholds list(thr_low, thr_high) used to build the C features (store what you used).
#' @param nFolds Number of cross-validation folds.
#' @param alpha Elastic-net mixing parameter.
#' @param seed Random seed.
#' @return A model usable by [predict_sysHIX()] / [compute_S_C()], plus `training` (sample scores).
#' @note The returned training sysHIX is computed on the samples the fusion model was fitted on.
#' @export
fit_sysHIX_model <- function(tab, levels = c("NT", "Inflammatory", "Malignant"), scHIX_thresholds,
                             nFolds = 5, alpha = 0.1, seed = 123) {
  .need(c("ordinalNet", "caret"), "training the fusion model")
  if (!length(levels) %in% 2:3) stop("levels must have length 2 or 3.")
  tab <- tab[tab$label %in% levels, ]
  y <- factor(tab$label, levels = levels, ordered = TRUE)
  num <- tab[, vapply(tab, is.numeric, logical(1)), drop = FALSE]
  na_cols <- names(which(colSums(is.na(num)) > 0))
  if (length(na_cols)) stop("Missing values in training features: ", paste(na_cols, collapse = ", "))
  if (any(table(y) < 5)) warning("A class has fewer than 5 samples: ", paste(names(table(y)), table(y), collapse = ", "), call. = FALSE)
  if (is.null(tab$series)) warning("No `series` column: out-of-fold scores are not cohort-grouped and may be optimistic.", call. = FALSE)
  cm <- grep("^M_CM\\d+$", names(tab), value = TRUE)
  x_cols <- c("E_seHIX", "C_Low_Frac", "C_Mid_Frac", "C_High_Frac", "C_mean", cm)
  nt <- y == levels[1]
  if (sum(nt) < 10) warning("Only ", sum(nt), " NT samples for standardisation.", call. = FALSE)

  if (length(unique(tab$n_epi_cells)) < 3)
    warning("Epithelial cell number is (nearly) constant across samples; cell-number adjustment skipped.", call. = FALSE)
  res_params <- list()
  for (cc in x_cols) {
    if (length(unique(tab$n_epi_cells)) < 3) {   # cell number (almost) constant: nothing to regress out
      res_params[[cc]] <- list(intercept = mean(tab[[cc]]), slope = 0)
    } else {
      fit <- stats::lm(tab[[cc]] ~ log1p(tab$n_epi_cells))
      res_params[[cc]] <- list(intercept = unname(stats::coef(fit)[1]), slope = unname(stats::coef(fit)[2]))
    }
    tab[[paste0(cc, "_res")]] <- tab[[cc]] - (res_params[[cc]]$intercept + res_params[[cc]]$slope * log1p(tab$n_epi_cells))
  }
  zs <- function(cols) { P <- list()
    for (f in cols) { P[[f]] <- list(mean = mean(tab[[f]][nt], na.rm = TRUE), sd = stats::sd(tab[[f]][nt], na.rm = TRUE) + 1e-8)
      tab[[paste0(f, "_z")]] <<- (tab[[f]] - P[[f]]$mean) / P[[f]]$sd }
    P }
  C_res <- paste0(c("C_Low_Frac", "C_Mid_Frac", "C_High_Frac", "C_mean"), "_res"); M_res <- paste0(cm, "_res")
  C_norm <- zs(C_res); M_norm <- zs(M_res)
  XC <- as.matrix(tab[, paste0(C_res, "_z")]); XM <- as.matrix(tab[, paste0(M_res, "_z")])
  set.seed(seed)
  lC <- .fit_layer(XC, y, nFolds, alpha); lM <- .fit_layer(XM, y, nFolds, alpha)
  folds <- .grouped_folds(y, tab$series, nFolds, seed)
  tab$S_C <- .oof_score(XC, y, lC$lambda, folds, alpha); tab$S_M <- .oof_score(XM, y, lM$lambda, folds, alpha)

  E_mu <- mean(tab$E_seHIX_res[nt], na.rm = TRUE); E_sd <- stats::sd(tab$E_seHIX_res[nt], na.rm = TRUE) + 1e-8
  tab$E_seHIX_z <- (tab$E_seHIX_res - E_mu) / E_sd
  SC_mu <- mean(tab$S_C[nt]); SC_sd <- stats::sd(tab$S_C[nt]) + 1e-8; tab$S_C_z <- (tab$S_C - SC_mu) / SC_sd
  SM_mu <- mean(tab$S_M[nt]); SM_sd <- stats::sd(tab$S_M[nt]) + 1e-8; tab$S_M_z <- (tab$S_M - SM_mu) / SM_sd
  fz <- c("E_seHIX_z", "S_C_z", "S_M_z")
  lF <- .fit_layer(as.matrix(tab[, fz]), y, nFolds, alpha)

  mapping <- stats::setNames(seq_along(levels) - 1, levels)
  model <- list(
    metadata = list(model_version = "user-trained", creation_date = Sys.Date(), R_version = R.version.string,
                    severity_levels = levels, severity_mapping = mapping, n_train = nrow(tab), n_classes = length(levels)),
    C_layer = list(model = lC$fit, beta = lC$beta, lambda = lC$lambda, norm_params = C_norm,
                   features = paste0(C_res, "_z"), features_raw = sub("_res$", "", C_res)),
    M_layer = list(model = lM$fit, beta = lM$beta, lambda = lM$lambda, norm_params = M_norm,
                   features = paste0(M_res, "_z"), features_raw = sub("_res$", "", M_res)),
    final_layer = list(model = lF$fit, beta = lF$beta, lambda = lF$lambda, features = fz),
    E_layer = list(mean = E_mu, sd = E_sd),
    layer_scores_norm = list(S_C = list(mean = SC_mu, sd = SC_sd), S_M = list(mean = SM_mu, sd = SM_sd)),
    scHIX_thresholds = c(scHIX_thresholds[c("thr_low", "thr_high")]),
    res_params = res_params)
  cum <- stats::predict(lF$fit, newx = as.matrix(tab[, fz]), type = "response", lambda = lF$lambda)
  tab$sysHIX <- as.numeric(.cum_to_class(cum, levels) %*% mapping)
  attr(model, "training") <- tab[, c("sample_id", "label", intersect("series", names(tab)), fz, "sysHIX")]
  class(model) <- c("sysHIX_model", "list")
  model
}
