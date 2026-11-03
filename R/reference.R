# Annotation references: download/cache, loading, validation and the marker-gene list used at the second level

.reference_registry <- list(
  HIXscore = list(file = "scRNA_ref_HIXscore.rds",
                  url = "https://comics.med.sustech.edu.cn/track/scRNA-ref/scRNA_ref_HIXscore.rds",
                  md5 = "39338840816b5d78d08fd89fe2a5e46c", size = 552876354)
)

.reference_cache_dir <- function() {
  d <- tools::R_user_dir("HIXscore", which = "cache")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

#' Marker genes of the second annotation level
#'
#' Named list (one element per major cell type) with the marker genes of the cell states of that type. They are
#' added to the genes used to transfer cell-state labels in [annotate_cells()].
#' @return Named list of character vectors.
#' @export
default_ref_features <- function() {
  readRDS(system.file("extdata", "ref_features_list.rds", package = "HIXscore", mustWork = TRUE))
}

#' Save an annotation reference
#' @param reference Object from [prepare_reference()].
#' @param path Output `.rds` file.
#' @return `path`, invisibly, with the md5 sum as attribute `md5`.
#' @export
save_reference <- function(reference, path) {
  stopifnot(inherits(reference, "HIXscore_reference"))
  saveRDS(reference, path)
  invisible(structure(path, md5 = unname(tools::md5sum(path))))
}

#' Download the annotation reference to the cache
#'
#' Downloads the HIXscore single-cell reference (a Seurat object of about 550 MB) into the user cache
#' (`tools::R_user_dir("HIXscore", "cache")`) and verifies its md5 sum. A cached file with the correct md5 sum is
#' reused. The md5 sum detects a damaged or changed download; it does not establish who produced the file.
#' @param force Download again even if a verified copy is cached.
#' @param timeout Download time limit in seconds.
#' @return Path of the cached file.
#' @export
download_reference <- function(force = FALSE, timeout = 3600) {
  r <- .reference_registry$HIXscore; dest <- file.path(.reference_cache_dir(), r$file)
  if (file.exists(dest) && !force && unname(tools::md5sum(dest)) == r$md5) return(dest)
  old <- options(timeout = max(timeout, getOption("timeout"))); on.exit(options(old))
  tmp <- paste0(dest, ".part"); on.exit(unlink(tmp), add = TRUE)
  utils::download.file(r$url, tmp, mode = "wb", quiet = FALSE)
  if (unname(tools::md5sum(tmp)) != r$md5) stop("md5 mismatch for the downloaded reference; file removed.", call. = FALSE)
  file.rename(tmp, dest)
  dest
}

#' Check that a reference covers the cell states of an M-layer model
#' @param reference A `HIXscore_reference`.
#' @param model M-layer model (default shipped model).
#' @return Invisibly, the model states missing from the reference (warns if any).
#' @export
check_reference <- function(reference, model = default_M_model()) {
  stopifnot(inherits(reference, "HIXscore_reference"))
  miss <- setdiff(rownames(model$W), reference$states)
  if (length(miss)) warning(length(miss), " of ", nrow(model$W), " M-layer cell states are absent from the reference: ",
                            paste(utils::head(miss, 5), collapse = ", "), if (length(miss) > 5) ", ...", call. = FALSE)
  invisible(miss)
}

# validate a raw reference Seurat object before it is prepared
.validate_reference_object <- function(ref, major_col, minor_col, reduction) {
  if (!inherits(ref, "Seurat")) stop("The reference must be a Seurat object, not ", class(ref)[1], ".", call. = FALSE)
  md <- ref@meta.data
  miss <- setdiff(c(major_col, minor_col), names(md))
  if (length(miss)) stop("The reference lacks metadata column(s): ", paste(miss, collapse = ", "), ".", call. = FALSE)
  if (anyNA(md[[major_col]]) || anyNA(md[[minor_col]])) stop("The reference has missing values in '", major_col, "' or '", minor_col, "'.", call. = FALSE)
  if (!"Epithelial" %in% md[[major_col]]) stop("The reference has no cells labelled 'Epithelial' in '", major_col, "'.", call. = FALSE)
  if (!"RNA" %in% names(ref@assays)) stop("The reference has no 'RNA' assay.", call. = FALSE)
  cnt <- ref@assays$RNA@counts
  if (length(cnt) == 0 || length(cnt@x) == 0) stop("The reference 'RNA' assay has no raw counts.", call. = FALSE)
  if (!reduction %in% names(ref@reductions)) stop("The reference has no '", reduction, "' reduction (available: ", paste(names(ref@reductions), collapse = ", "), ").", call. = FALSE)
  L <- Seurat::Loadings(ref, reduction)
  if (nrow(L) == 0) stop("The '", reduction, "' reduction of the reference has no feature loadings.", call. = FALSE)
  if (nrow(Seurat::Embeddings(ref, reduction)) != ncol(ref)) stop("The '", reduction, "' embeddings do not cover all cells of the reference.", call. = FALSE)
  invisible(TRUE)
}

#' Load the annotation reference
#'
#' Reads a reference, checks it and prepares it for [annotate_cells()]. Without `path` the HIXscore reference is
#' downloaded to the cache ([download_reference()]). `path` can point to the downloaded Seurat object, to your own
#' reference (see [prepare_reference()] for what it has to contain), or to a prepared reference saved with
#' [save_reference()], which is returned as it is.
#'
#' Preparing the reference splits it by major cell type and recomputes a PCA within each type, so it needs a
#' few GB of memory for the 50,000-cell HIXscore reference. Reading an `.rds` file runs R code stored in the file:
#' only load files from sources you trust.
#'
#' @param path Optional path to an `.rds` file.
#' @param ref_features Marker genes of the second level, default [default_ref_features()] (named by major type).
#' @param reduction Reduction used for label transfer (default `"pca"`).
#' @param major_col,minor_col Metadata columns with the major cell type and the cell state.
#' @param min_cells_per_type Minimum reference cells of a cell state.
#' @return A `HIXscore_reference`.
#' @export
load_reference <- function(path = NULL, ref_features = default_ref_features(), reduction = "pca",
                           major_col = "celltype_major", minor_col = "cell_type2", min_cells_per_type = 50) {
  .check_seurat()
  if (is.null(path)) path <- download_reference()
  if (!file.exists(path)) stop("Reference file not found: ", path, call. = FALSE)
  ref <- readRDS(path)
  if (inherits(ref, "HIXscore_reference")) { check_reference(ref); return(ref) }
  .validate_reference_object(ref, major_col, minor_col, reduction)
  ref <- prepare_reference(ref, ref_features = ref_features, major_col = major_col, minor_col = minor_col,
                           reduction = reduction, min_cells_per_type = min_cells_per_type)
  check_reference(ref)
  ref
}
