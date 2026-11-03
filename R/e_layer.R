utils::globalVariables(c("sample_id", "tissue", "label", "value", "group"))

# E layer: seHIX (cosine distance of MP activity + heterogeneity profile to the NT baseline)

.mp_cols <- function(md, mp_pattern = "^MP_\\d+_", exclude_na_mp = TRUE) {
  cols <- grep(mp_pattern, colnames(md), value = TRUE)
  if (exclude_na_mp) cols <- cols[!grepl("_NA$|_Ribosome", cols)]
  if (length(cols) == 0) stop("No MP columns matching '", mp_pattern, "' found in the metadata.")
  cols
}

# sample-level features: mean activity and CV (heterogeneity) of each MP
.e_sample_features <- function(md, mp_cols, sample_col, tissue_col, min_cells) {
  stopifnot(sample_col %in% colnames(md), tissue_col %in% colnames(md))
  md <- md[, c(sample_col, tissue_col, mp_cols)]
  names(md)[1:2] <- c("sample_id", "tissue")
  act <- dplyr::summarise(dplyr::group_by(md, sample_id, tissue),
                          dplyr::across(dplyr::all_of(mp_cols), ~ mean(.x, na.rm = TRUE)),
                          n_cells = dplyr::n(), .groups = "drop")
  cv <- dplyr::summarise(dplyr::group_by(md, sample_id, tissue),
                         dplyr::across(dplyr::all_of(mp_cols),
                                       ~ stats::sd(.x, na.rm = TRUE) / (mean(.x, na.rm = TRUE) + 1e-10),
                                       .names = "{.col}_CV"), .groups = "drop")
  out <- dplyr::inner_join(act, cv, by = c("sample_id", "tissue"))
  out[out$n_cells >= min_cells, , drop = FALSE]
}

.scale_with <- function(x, mu, sdv) (x - mu) / (sdv + 1e-10)

#' Fit an E-layer (seHIX) NT baseline
#'
#' Stores the normal-tissue (NT) reference vector together with the per-feature
#' mean/SD used to standardise sample features, so that new samples can be
#' scored independently of the other samples in their own dataset.
#'
#' @param meta Cell-level metadata (data.frame) with MP score columns, e.g. `seurat@meta.data`.
#' @param sample_col,tissue_col Metadata columns with sample id and tissue label.
#' @param nt_label Label of normal tissue in `tissue_col`.
#' @param min_cells Minimum epithelial cells per sample (default 30).
#' @param min_nt_samples Minimum number of NT samples (passing `min_cells`); fewer is an error.
#'   Between `min_nt_samples` and 9 NT samples a warning is given (the baseline and the feature
#'   standardisation are estimated from few samples).
#' @param mp_pattern Regular expression selecting MP columns.
#' @param exclude_na_mp Drop MPs labelled NA or Ribosome.
#' @return An object of class `sysHIX_E_baseline`.
#' @export
fit_E_baseline <- function(meta, sample_col = "sample_id", tissue_col = "tissue",
                           nt_label = "NT", mp_pattern = "^MP_\\d+_",
                           exclude_na_mp = TRUE, min_cells = 30, min_nt_samples = 3) {
  mp_cols <- .mp_cols(meta, mp_pattern, exclude_na_mp)
  f <- .e_sample_features(meta, mp_cols, sample_col, tissue_col, min_cells)
  feats <- c(mp_cols, paste0(mp_cols, "_CV"))
  n_nt <- sum(f$tissue == nt_label)
  if (n_nt < min_nt_samples)
    stop("Only ", n_nt, " NT sample(s) (tissue == '", nt_label, "') with >= ", min_cells, " epithelial cells; at least ",
         min_nt_samples, " are needed to fit an E-layer baseline.", call. = FALSE)
  if (n_nt < 10)
    warning("The E-layer baseline is fitted on only ", n_nt, " NT samples (< 10); E-layer scores are unstable. ",
            "Use the shipped oral baseline (mode = 'pretrained') or more NT samples if possible.", call. = FALSE)
  mu  <- vapply(feats, function(k) mean(f[[k]], na.rm = TRUE), numeric(1))
  sdv <- vapply(feats, function(k) stats::sd(f[[k]], na.rm = TRUE), numeric(1))
  scaled <- mapply(.scale_with, f[feats], mu, sdv, SIMPLIFY = FALSE)
  nt <- f$tissue == nt_label
  nt_vec <- vapply(scaled, function(v) mean(v[nt], na.rm = TRUE), numeric(1))
  structure(list(mp_cols = mp_cols, features = feats, mean = mu, sd = sdv,
                 nt_baseline = nt_vec, n_nt = sum(nt), min_cells = min_cells,
                 creation_date = Sys.Date()), class = "sysHIX_E_baseline")
}

# Accept the legacy rds layout written by save_E_seHIX_baseline()
.as_E_baseline <- function(b) {
  if (inherits(b, "sysHIX_E_baseline")) return(b)
  if (is.character(b) && length(b) == 1) b <- readRDS(b)
  if (inherits(b, "sysHIX_E_baseline")) return(b)
  stopifnot(is.list(b), !is.null(b$nt_baseline), !is.null(b$mp_cols))
  mp <- b$mp_cols
  feats <- c(mp, paste0(mp, "_CV"))
  np <- cbind(as.data.frame(b$activity_norm_params), as.data.frame(b$cv_norm_params))
  mu  <- vapply(feats, function(k) np[[paste0(k, "_mean")]], numeric(1))
  sdv <- vapply(feats, function(k) np[[paste0(k, "_sd")]], numeric(1))
  structure(list(mp_cols = mp, features = feats, mean = mu, sd = sdv,
                 nt_baseline = as.numeric(b$nt_baseline), n_nt = NA_integer_,
                 min_cells = 10, creation_date = b$creation_date), class = "sysHIX_E_baseline")
}

.cosine_dist <- function(x, y) {
  if (anyNA(x) || anyNA(y)) return(NA_real_)
  nx <- sqrt(sum(x^2)); ny <- sqrt(sum(y^2))
  if (nx < 1e-10 || ny < 1e-10) return(NA_real_)
  1 - sum(x * y) / (nx * ny)
}

#' Compute the E-layer score (seHIX)
#'
#' @param meta Cell-level metadata with MP score columns (epithelial cells only).
#' @param baseline `sysHIX_E_baseline`, a path to a baseline rds, or `NULL`.
#'   If `NULL`, the baseline is fitted on the NT samples of `meta` (organ-specific mode).
#' @param legacy_scaling If `TRUE` and an external `baseline` is given, standardise
#'   features using the mean/SD of the input samples. Default `FALSE`: use the frozen
#'   baseline mean/SD, so a sample's score does not depend on the other samples in the input.
#' @param min_cells Minimum epithelial cells per sample (default 30).
#' @param sample_col Metadata column with sample id.
#' @param tissue_col Metadata column with tissue label.
#' @param nt_label Label of normal tissue in the tissue column.
#' @param mp_pattern Regular expression selecting MP columns.
#' @param exclude_na_mp Drop MPs labelled NA or Ribosome.
#' @return data.frame with `sample_id`, `tissue`, `n_cells`, `E_seHIX`.
#' @export
compute_E_seHIX <- function(meta, baseline = NULL, sample_col = "sample_id",
                            tissue_col = "tissue", nt_label = "NT",
                            mp_pattern = "^MP_\\d+_", exclude_na_mp = TRUE,
                            legacy_scaling = FALSE, min_cells = 30) {
  if (is.null(baseline)) {
    baseline <- fit_E_baseline(meta, sample_col, tissue_col, nt_label, mp_pattern,
                               exclude_na_mp, min_cells)
  } else {
    baseline <- .as_E_baseline(baseline)
  }
  missing_mp <- setdiff(baseline$mp_cols, colnames(meta))
  if (length(missing_mp)) stop("Metadata lacks MP columns required by the baseline: ",
                               paste(missing_mp, collapse = ", "))
  f <- .e_sample_features(meta, baseline$mp_cols, sample_col, tissue_col, min_cells)
  if (nrow(f) == 0) stop("No sample has at least ", min_cells, " epithelial cells.")
  feats <- baseline$features
  if (legacy_scaling) {
    if (nrow(f) < 5) warning("legacy_scaling with fewer than 5 samples gives unstable scores.")
    mu  <- vapply(feats, function(k) mean(f[[k]], na.rm = TRUE), numeric(1))
    sdv <- vapply(feats, function(k) stats::sd(f[[k]], na.rm = TRUE), numeric(1))
  } else {
    mu <- baseline$mean; sdv <- baseline$sd
  }
  X <- vapply(feats, function(k) .scale_with(f[[k]], mu[[k]], sdv[[k]]), numeric(nrow(f)))
  X <- matrix(X, nrow = nrow(f), dimnames = list(NULL, feats))
  f$E_seHIX <- apply(X, 1, .cosine_dist, y = baseline$nt_baseline)
  out <- f[, c("sample_id", "tissue", "n_cells", "E_seHIX")]
  if (anyNA(out$E_seHIX)) {
    warning(sum(is.na(out$E_seHIX)), " sample(s) with zero or missing feature vector removed.")
    out <- out[!is.na(out$E_seHIX), ]
  }
  attr(out, "scaling") <- if (legacy_scaling) "input-dataset" else "frozen-baseline"
  attr(out, "n_nt_baseline") <- baseline$n_nt
  as.data.frame(out)
}
