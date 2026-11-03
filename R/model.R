# Fusion model: C-layer (S_C), M-layer (S_M) and final ordinal fusion (sysHIX)

#' Pretrained oral sysHIX model
#' @return Model list (C layer, M layer, ordinal final layer, normalisation parameters).
#' @export
default_sysHIX_model <- function() {
  readRDS(system.file("extdata", "sysHIX_model_oral.rds", package = "HIXscore", mustWork = TRUE))
}

# Cell-number effect removal using regression coefficients fitted on the training set only
.apply_cellcount_residual <- function(feature_name, new_data, model) {
  rp <- model$res_params[[feature_name]]
  if (is.null(rp))
    stop("Model lacks res_params for '", feature_name, "'; it was probably saved by an old version.")
  new_data[[feature_name]] - (rp$intercept + rp$slope * log1p(new_data$n_epi_cells))
}

.sample_C_features <- function(md, scHIX_col, sample_col, thresholds) {
  stopifnot(scHIX_col %in% colnames(md), sample_col %in% colnames(md))
  x <- md[[scHIX_col]]; s <- md[[sample_col]]
  do.call(rbind, lapply(split(seq_along(x), s), function(i) {
    v <- x[i]
    data.frame(n_epi_cells = length(i),
               C_Low_Frac  = mean(v <= thresholds$thr_low, na.rm = TRUE),
               C_Mid_Frac  = mean(v > thresholds$thr_low & v <= thresholds$thr_high, na.rm = TRUE),
               C_High_Frac = mean(v > thresholds$thr_high, na.rm = TRUE),
               C_mean      = stats::median(v, na.rm = TRUE))
  })) -> out
  out <- cbind(sample_id = rownames(out), out, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

.layer_z <- function(df, layer, model) {
  raw <- gsub("_res$", "", names(model[[layer]]$norm_params))
  for (f in raw) {
    fr <- paste0(f, "_res")
    df[[fr]] <- .apply_cellcount_residual(f, df, model)
    p <- model[[layer]]$norm_params[[fr]]
    df[[paste0(fr, "_z")]] <- (df[[fr]] - p$mean) / p$sd
  }
  df
}

#' Compute the C-layer score S_C
#'
#' @param x Seurat object or cell-level metadata of epithelial cells with a scHIX column.
#' @param model sysHIX model, default [default_sysHIX_model()].
#' @param min_cells Minimum epithelial cells per sample (default 30).
#' @param scHIX_col Metadata column with cell-level scHIX.
#' @param sample_col Metadata column with sample id.
#' @return data.frame with `sample_id`, `n_epi_cells`, C features, `S_C`, `S_C_z`.
#' @export
compute_S_C <- function(x, model = default_sysHIX_model(), scHIX_col = "scHIX_cosine",
                        sample_col = "sample_id", min_cells = 30) {
  f <- .sample_C_features(.get_meta(x), scHIX_col, sample_col, model$scHIX_thresholds)
  f <- f[f$n_epi_cells >= min_cells, , drop = FALSE]
  if (nrow(f) == 0) stop("No sample has at least ", min_cells, " epithelial cells.")
  .S_C_from_features(f, model)
}

# sample-level C features (n_epi_cells, C_*) -> S_C, S_C_z with a model's frozen parameters
.S_C_from_features <- function(f, model) {
  f <- .layer_z(f, "C_layer", model)
  X <- as.matrix(f[, model$C_layer$features, drop = FALSE])
  f$S_C <- -as.numeric(X %*% model$C_layer$beta)
  f$S_C_z <- (f$S_C - model$layer_scores_norm$S_C$mean) / model$layer_scores_norm$S_C$sd
  f
}

.cum_to_class <- function(cum, levels) {
  K <- length(levels); n <- nrow(cum)
  p <- matrix(0, n, K, dimnames = list(NULL, levels))
  p[, 1] <- cum[, 1]
  if (K > 2) for (k in 2:(K - 1)) p[, k] <- cum[, k] - cum[, k - 1]
  p[, K] <- 1 - cum[, K - 1]
  p[p < 0] <- 0
  p / rowSums(p)
}

#' E-layer standardised score
#'
#' Removes the epithelial cell-number effect with the training-set regression and standardises
#' with the frozen training parameters. Uses only the E-layer samples.
#' @param E Output of [compute_E_seHIX()] (needs `E_seHIX` and `n_cells`).
#' @param model sysHIX model, default [default_sysHIX_model()].
#' @return `E` with columns `E_seHIX_res` and `E_seHIX_z`.
#' @export
compute_E_z <- function(E, model = default_sysHIX_model()) {
  E_in <- E
  d <- data.frame(E_seHIX = E$E_seHIX, n_epi_cells = E$n_cells)
  E$E_seHIX_res <- .apply_cellcount_residual("E_seHIX", d, model)
  E$E_seHIX_z <- (E$E_seHIX_res - model$E_layer$mean) / model$E_layer$sd
  attr(E, "n_nt_baseline") <- attr(E_in, "n_nt_baseline")
  E
}

#' M-layer score S_M
#'
#' Combines the CM scores into S_M with the frozen model. The model removes the effect of the
#' epithelial cell number of the same sample from every CM, so `n_epi` (all epithelial cells of
#' each sample, **before** any epithelial cell-number filter) is needed. Samples without epithelial
#' cells are kept (`n_epi` = 0). Samples with at most `min_epi_training` epithelial cells are flagged in
#' `M_epi_below_training_min`: their epithelial cell-number adjustment lies below the range of the
#' training samples.
#'
#' @param M Output of [compute_M_CM_scores()].
#' @param n_epi Named numeric vector (names = sample ids) of epithelial cells per sample.
#' @param model sysHIX model, default [default_sysHIX_model()].
#' @param min_epi_training Minimum epithelial cells of the training samples, used only for the flag.
#' @param min_epi Samples with fewer epithelial cells are dropped; default 0 keeps all.
#' @return `M` with `n_epi_cells`, `S_M`, `S_M_z` and the logical `M_epi_below_training_min`.
#' @export
compute_S_M <- function(M, n_epi, model = default_sysHIX_model(), min_epi_training = 10, min_epi = 0) {
  Mraw <- gsub("_res$", "", names(model$M_layer$norm_params))
  if (!all(Mraw %in% names(M))) stop("M scores lack columns required by the model: ", paste(setdiff(Mraw, names(M)), collapse = ", "))
  d <- M
  d$n_epi_cells <- as.numeric(n_epi[as.character(d$sample_id)]); d$n_epi_cells[is.na(d$n_epi_cells)] <- 0
  d <- d[d$n_epi_cells >= min_epi, , drop = FALSE]
  if (nrow(d) == 0) stop("No sample left after the epithelial cell-number filter for the M layer.")
  d <- .layer_z(d, "M_layer", model)
  d$S_M <- -as.numeric(as.matrix(d[, model$M_layer$features, drop = FALSE]) %*% model$M_layer$beta)
  d$S_M_z <- (d$S_M - model$layer_scores_norm$S_M$mean) / model$layer_scores_norm$S_M$sd
  d$M_epi_below_training_min <- d$n_epi_cells <= min_epi_training
  d[, !grepl("_res(_z)?$", names(d)), drop = FALSE]
}

#' Fuse the three layers into sysHIX
#'
#' Each layer is scored independently on its own samples; **only here** the samples are intersected:
#' sysHIX is computed for samples that have all of `E_seHIX_z`, `S_C_z` and `S_M_z`.
#'
#' @param E_scores Output of [compute_E_z()] (or of [compute_E_seHIX()], then `E_seHIX_z` is added).
#' @param C_scores Output of [compute_S_C()].
#' @param M_scores Output of [compute_S_M()] (or of [compute_M_CM_scores()], then `S_M_z` is added
#'   using the epithelial cell counts in `C_scores`).
#' @param model sysHIX model, default [default_sysHIX_model()].
#' @return data.frame, one row per E, C and M sample: `E_seHIX_z`, `S_C_z`, `S_M_z`, `sysHIX`, class probabilities.
#' @export
predict_sysHIX <- function(E_scores, C_scores, M_scores, model = default_sysHIX_model()) {
  .need("ordinalNet", "predict_sysHIX()")
  if (!"E_seHIX_z" %in% names(E_scores)) E_scores <- compute_E_z(E_scores, model)
  if (!"S_M_z" %in% names(M_scores))
    M_scores <- compute_S_M(M_scores, stats::setNames(C_scores$n_epi_cells, C_scores$sample_id), model)
  d <- merge(E_scores[, c("sample_id", "E_seHIX_z")], C_scores[, c("sample_id", "S_C_z")], by = "sample_id")
  d <- merge(d, M_scores[, c("sample_id", "S_M_z")], by = "sample_id")
  if (nrow(d) == 0) stop("No sample has E, C and M scores.")
  X <- as.matrix(d[, model$final_layer$features])
  cum <- stats::predict(model$final_layer$model, newx = X, type = "response", lambda = model$final_layer$lambda)
  prob <- .cum_to_class(cum, model$metadata$severity_levels)
  d$sysHIX <- as.numeric(prob %*% model$metadata$severity_mapping)
  colnames(prob) <- paste0("prob_", colnames(prob))
  cbind(d, prob)
}
