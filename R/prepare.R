# From a raw (annotated or not) Seurat object to the two inputs of hixscore()

#' Prepare input for [hixscore()]
#'
#' Genes expressed in fewer than `min_cells_gene` cells are dropped, all cells are log-normalised
#' together, and the cells are split by major cell type into epithelial cells (`epi`, E and C layers)
#' and non-epithelial cells (`tme`, M layer).
#'
#' Cells must carry the labels `predicted.major` (with `"Epithelial"` for epithelial cells) and
#' `predicted.minor` (cell state, the vocabulary of the M-layer model). Either supply them, or set
#' `annotate = TRUE` and a `reference` to compute them with [annotate_cells()].
#'
#' @param obj Seurat object with raw counts.
#' @param annotate Annotate cells with [annotate_cells()] first.
#' @param reference A `HIXscore_reference`; required when `annotate = TRUE`.
#' @param sample_col Metadata column with the sample id.
#' @param major_col,minor_col Metadata columns with major type and cell state (default: the annotation output).
#' @param epithelial_label Value of `major_col` marking epithelial cells.
#' @param min_cells_gene Genes detected in fewer cells are removed before normalisation.
#' @param ... Passed to [annotate_cells()].
#' @return list(epi, tme, summary): `epi` Seurat object (log-normalised); `tme` data.frame with `sample_id`
#'   and `cell_state`; `summary` the annotation summary, including the fraction of non-epithelial cells whose
#'   state is in the M-layer model.
#' @export
prepare_input <- function(obj, annotate = FALSE, reference = NULL, sample_col = "sample_id",
                          major_col = "predicted.major", minor_col = "predicted.minor",
                          epithelial_label = "Epithelial", min_cells_gene = 30, ...) {
  .check_seurat()
  if (annotate) {
    if (is.null(reference)) stop("`annotate = TRUE` needs a `reference` (see load_reference() / prepare_reference()).", call. = FALSE)
    obj <- annotate_cells(obj, reference, sample_col = sample_col, ...)
    major_col <- "predicted.major"; minor_col <- "predicted.minor"
  }
  md <- obj@meta.data
  miss <- setdiff(c(sample_col, major_col, minor_col), names(md))
  if (length(miss)) stop("Missing metadata column(s): ", paste(miss, collapse = ", "),
                         ". Annotate the cells first (annotate = TRUE) or supply them.", call. = FALSE)
  cnt <- Seurat::GetAssayData(obj, slot = "counts")
  keep <- rownames(cnt)[Matrix::rowSums(cnt > 0) >= min_cells_gene]
  if (length(keep) < 1000) warning("Only ", length(keep), " genes pass the expression filter; MP/module scores may be unreliable.", call. = FALSE)
  obj <- Seurat::NormalizeData(subset(obj, features = keep), verbose = FALSE)
  is_epi <- !is.na(md[[major_col]]) & md[[major_col]] == epithelial_label
  if (!any(is_epi)) stop("No cells labelled '", epithelial_label, "' in '", major_col, "'.", call. = FALSE)
  non_epi <- !is.na(md[[major_col]]) & !is_epi
  tme <- data.frame(sample_id = as.character(md[[sample_col]][non_epi]), cell_state = as.character(md[[minor_col]][non_epi]),
                    stringsAsFactors = FALSE, row.names = rownames(md)[non_epi])
  sm <- annotation_summary(data.frame(predicted.major = md[[major_col]], predicted.minor = md[[minor_col]]))
  if (is.finite(sm$in_model_state_fraction) && sm$in_model_state_fraction < 0.5)
    warning("Only ", round(100 * sm$in_model_state_fraction), "% of non-epithelial cells have a cell state used by the M-layer model; ",
            "check that the labels follow the model vocabulary (rownames(default_M_model()$W)).", call. = FALSE)
  list(epi = obj[, rownames(md)[is_epi]], tme = tme, summary = sm)
}
