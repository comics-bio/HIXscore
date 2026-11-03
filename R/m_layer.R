# M layer: tumour-microenvironment ecotype (CM) scores from immune/stromal cell-state composition

.get_meta <- function(x) if (inherits(x, "Seurat")) x@meta.data else as.data.frame(x)

#' Pretrained oral M-layer model (NMF with k = 4 cell-state modules)
#'
#' A list with the NMF basis `W` (cell states x 4 CMs) and the per-state
#' minimum/maximum of the training composition (proportions), used to scale new samples.
#' @export
default_M_model <- function() {
  readRDS(system.file("extdata", "M_model_oral_K4.rds", package = "HIXscore", mustWork = TRUE))
}

#' Compute M-layer CM scores
#'
#' For each sample, the proportion of every immune/stromal cell state is computed
#' (denominator: cells assigned to the model's states), scaled to 0 to 1 with the
#' **frozen training** per-state min/max, and non-negative least squares gives the
#' CM scores (`V ~ W H`, W fixed), rescaled to sum to 1 per sample. A sample's score therefore does not depend on the
#' other samples in the input.
#'
#' @param x Seurat object or cell-level metadata of non-epithelial cells.
#' @param model M model from [default_M_model()] (or built by `build_M_ecotypes()`).
#' @param sample_col,state_col Metadata columns with sample id and cell-state label.
#' @param min_cells Minimum immune+stromal cells per sample (default 200).
#' @return data.frame with `sample_id`, `M_CM1`...`M_CMk`, `n_cells`, `frac_in_model_states`.
#' @export
compute_M_CM_scores <- function(x, model = default_M_model(), sample_col = "sample_id",
                                state_col = "cell_state", min_cells = 200) {
  md <- .get_meta(x)
  stopifnot(sample_col %in% colnames(md), state_col %in% colnames(md))
  md <- md[!is.na(md[[state_col]]), c(sample_col, state_col)]
  names(md) <- c("sample_id", "state")
  md$state <- sub("^Fibroblasts_", "Fibro_", as.character(md$state))
  W <- model$W; states <- rownames(W); k <- ncol(W)
  tot <- table(md$sample_id)
  keep <- names(tot)[tot >= min_cells]
  if (length(keep) == 0) stop("No sample has at least ", min_cells, " immune/stromal cells.")
  in_model <- md$state %in% states
  if (!any(in_model)) stop("None of the model's cell states were found in '", state_col, "'. ",
                           "Cell-state labels must use the model vocabulary (see rownames(model$W)).")
  md <- md[md$sample_id %in% keep, ]
  in_model <- md$state %in% states
  cnt <- table(factor(md$state[in_model], levels = states), md$sample_id[in_model])
  cnt <- cnt[, keep, drop = FALSE]
  prop <- sweep(unclass(cnt), 2, pmax(colSums(cnt), 1), "/")
  rng <- pmax(model$comp_max - model$comp_min, 1e-12)
  V <- (prop - model$comp_min[states]) / rng[states]
  H <- vapply(seq_len(ncol(V)), function(i) {
    h <- nnls::nnls(W, V[, i])$x
    if (sum(h) > 0) h / sum(h) else rep(NA_real_, k)   # CM scores are compositional (sum to 1)
  }, numeric(k))
  H <- matrix(H, nrow = k)
  out <- data.frame(sample_id = colnames(V), t(H), check.names = FALSE, row.names = NULL)
  names(out)[-1] <- paste0("M_CM", seq_len(k))
  bad <- is.na(out$M_CM1)
  if (any(bad)) { warning(sum(bad), " sample(s) with an all-zero NNLS solution removed."); }
  out$n_cells <- as.integer(tot[out$sample_id])
  out$frac_in_model_states <- as.numeric(colSums(cnt) / tot[out$sample_id])
  out <- out[!bad, , drop = FALSE]
  if (any(out$frac_in_model_states < 0.5))
    warning(sum(out$frac_in_model_states < 0.5), " sample(s) have <50% of cells in the model's cell states; ",
            "M-layer scores may be unreliable.", call. = FALSE)
  out
}
