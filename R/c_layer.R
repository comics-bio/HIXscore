# C layer: scHIX (cell-level cosine distance of BDM-module scores to the NT baseline)

#' Default BDM (baseline-deviation module) gene programs shipped with the package
#' @return Named list of 5 character vectors (Module_1 ... Module_5).
#' @export
default_BDM_programs <- function() {
  readRDS(system.file("extdata", "BDM_programs.rds", package = "HIXscore", mustWork = TRUE))
}

#' Validate user-supplied BDM gene programs
#'
#' Checks that `programs` is a named list of non-empty, unique character gene vectors with unique
#' names, and (optionally) reports the overlap with the genes of your data.
#' @param programs Named list of gene sets.
#' @param min_genes Minimum genes per module.
#' @param genes Optional character vector, e.g. `rownames(obj)`.
#' @return `programs` (invisibly), with duplicated genes within a module removed.
#' @export
validate_BDM_programs <- function(programs, genes = NULL, min_genes = 5) {
  if (!is.list(programs) || length(programs) < 2) stop("`programs` must be a list of at least 2 gene sets.")
  if (is.null(names(programs)) || anyNA(names(programs)) || any(names(programs) == "") || anyDuplicated(names(programs)))
    stop("`programs` must have unique, non-empty names (e.g. Module_1, Module_2, ...).")
  if (!all(vapply(programs, is.character, logical(1)))) stop("Each gene set must be a character vector of gene symbols.")
  programs <- lapply(programs, function(g) unique(g[!is.na(g) & g != ""]))
  small <- vapply(programs, length, integer(1)) < min_genes
  if (any(small)) stop("Gene sets with fewer than ", min_genes, " genes: ", paste(names(programs)[small], collapse = ", "))
  if (!is.null(genes)) .check_program_overlap(programs, genes)
  invisible(programs)
}

.check_program_overlap <- function(programs, genes, min_frac = 0.7) {
  frac <- vapply(programs, function(g) mean(g %in% genes), numeric(1))
  if (any(frac < min_frac))
    warning("Low gene overlap with the BDM programs (", paste0(names(frac), "=", round(frac, 2), collapse = ", "),
            "). Scores may be unreliable.", call. = FALSE)
  invisible(frac)
}

.module_scores <- function(obj, programs, assay) {
  .check_program_overlap(programs, rownames(obj))
  tmp <- "sysHIXtmpmod"
  obj <- Seurat::AddModuleScore(obj, features = programs, name = tmp, assay = assay)
  cols <- paste0(tmp, seq_along(programs))
  m <- as.matrix(obj@meta.data[, cols, drop = FALSE])
  obj@meta.data[, cols] <- NULL
  colnames(m) <- names(programs)
  list(obj = obj, scores = m)
}

#' Fit a C-layer (scHIX) NT baseline
#'
#' @param obj Seurat object of epithelial cells (normalised).
#' @param programs Named list of gene programs; default [default_BDM_programs()].
#' @param tissue_col,nt_label Metadata column / label identifying NT cells.
#' @param min_nt_cells Minimum number of NT cells.
#' @param assay Seurat assay.
#' @return Object of class `sysHIX_C_baseline`.
#' @export
fit_C_baseline <- function(obj, programs = default_BDM_programs(), tissue_col = "tissue",
                           nt_label = "NT", assay = "RNA", min_nt_cells = 10) {
  programs <- validate_BDM_programs(programs)
  ms <- .module_scores(obj, programs, assay)
  nt <- ms$obj@meta.data[[tissue_col]] == nt_label
  if (sum(nt, na.rm = TRUE) < min_nt_cells)
    stop("Too few NT cells: ", sum(nt, na.rm = TRUE), " < ", min_nt_cells)
  structure(list(nt_baseline = colMeans(ms$scores[which(nt), , drop = FALSE]),
                 programs = programs, n_nt_cells = sum(nt), creation_date = Sys.Date()),
            class = "sysHIX_C_baseline")
}

.as_C_baseline <- function(b) {
  if (inherits(b, "sysHIX_C_baseline")) return(b)
  if (is.character(b) && length(b) == 1) b <- readRDS(b)
  if (inherits(b, "sysHIX_C_baseline")) return(b)
  stopifnot(is.list(b), !is.null(b$nt_baseline))
  structure(list(nt_baseline = b$nt_baseline, programs = b$baseline_programs,
                 n_nt_cells = b$n_nt_cells, creation_date = b$creation_date),
            class = "sysHIX_C_baseline")
}

.dist_fun <- function(method) switch(method,
  cosine = function(x, y) 1 - sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2))),
  euclidean = function(x, y) sqrt(sum((x - y)^2)),
  manhattan = function(x, y) sum(abs(x - y)),
  stop("Unknown distance method: ", method))

#' Compute cell-level scHIX
#'
#' Scores each cell on the BDM programs (`Seurat::AddModuleScore`) and computes its
#' distance to the NT baseline. `AddModuleScore` selects control genes using expression bins of the
#' input object.
#'
#' @param obj Seurat object of epithelial cells.
#' @param baseline `sysHIX_C_baseline`, a path to an rds, or `NULL` (fit on NT cells of `obj`).
#' @param programs Gene programs; taken from `baseline` when one is supplied.
#' @param tissue_col Metadata column with tissue label.
#' @param nt_label Label of normal tissue in the tissue column.
#' @param assay Seurat assay.
#' @param distance_method Distance to the NT baseline.
#' @param min_nt_cells Minimum NT cells.
#' @return `obj` with metadata column `scHIX_<distance_method>`.
#' @export
compute_scHIX <- function(obj, baseline = NULL, programs = NULL, tissue_col = "tissue",
                          nt_label = "NT", assay = "RNA", distance_method = "cosine",
                          min_nt_cells = 10) {
  if (is.null(baseline)) {
    if (is.null(programs)) programs <- default_BDM_programs()
    baseline <- fit_C_baseline(obj, programs, tissue_col, nt_label, assay, min_nt_cells)
  } else {
    baseline <- .as_C_baseline(baseline)
    if (!is.null(programs)) baseline$programs <- validate_BDM_programs(programs)
  }
  if (length(baseline$nt_baseline) != length(baseline$programs))
    stop("Baseline dimension does not match the number of programs.")
  ms <- .module_scores(obj, baseline$programs, assay)
  d <- .dist_fun(distance_method)
  col <- paste0("scHIX_", distance_method)
  ms$obj@meta.data[[col]] <- apply(ms$scores, 1, d, y = baseline$nt_baseline)
  ms$obj
}

#' scHIX tertile thresholds from NT cells
#'
#' Used when the NT baseline is refitted (mode `custom_baseline`): the low/high cut-offs
#' of the C-layer features are the 1/3 and 2/3 quantiles of the cell-level scHIX of NT cells.
#' @param x Seurat object or cell-level metadata.
#' @param scHIX_col Metadata column with cell-level scHIX.
#' @param tissue_col Metadata column with tissue label.
#' @param nt_label Label of normal tissue in the tissue column.
#' @param min_nt_cells Minimum NT cells.
#' @return list(thr_low, thr_high, n_nt_cells)
#' @export
nt_scHIX_thresholds <- function(x, scHIX_col = "scHIX_cosine", tissue_col = "tissue",
                                nt_label = "NT", min_nt_cells = 10) {
  md <- .get_meta(x)
  v <- md[[scHIX_col]][md[[tissue_col]] == nt_label]
  v <- v[!is.na(v)]
  if (length(v) < min_nt_cells) stop("Too few NT cells: ", length(v), " < ", min_nt_cells)
  q <- stats::quantile(v, c(1/3, 2/3), na.rm = TRUE)
  list(thr_low = unname(q[1]), thr_high = unname(q[2]), n_nt_cells = length(v))
}
