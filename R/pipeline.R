# Entry points: MP scoring, input checks and the one-call pipeline

#' Epithelial meta-program (MP) gene lists shipped with the package
#' @export
default_MP_list <- function() {
  readRDS(system.file("extdata", "MP_list.rds", package = "HIXscore", mustWork = TRUE))
}

#' Validate user-supplied E-layer gene sets (MPs)
#'
#' Checks that `mp_list` is a named list of character gene vectors with unique, non-empty, syntactically valid names
#' and at least `min_genes` genes per set, and removes duplicated or missing genes within a set.
#' @param mp_list Named list of gene sets.
#' @param min_genes Minimum genes per set.
#' @param genes Optional character vector, e.g. `rownames(obj)`; sets with low overlap give a warning.
#' @return `mp_list` (invisibly) with duplicated genes within a set removed.
#' @export
validate_MP_list <- function(mp_list, genes = NULL, min_genes = 5) {
  if (!is.list(mp_list) || length(mp_list) < 2) stop("`mp_list` must be a list of at least 2 gene sets.", call. = FALSE)
  nm <- names(mp_list)
  if (is.null(nm) || anyNA(nm) || any(nm == "") || anyDuplicated(nm))
    stop("`mp_list` must have unique, non-empty names.", call. = FALSE)
  if (any(nm != make.names(nm))) stop("`mp_list` names must be syntactically valid (letters, digits, '.', '_'; not starting with a digit): ", paste(nm[nm != make.names(nm)], collapse = ", "), call. = FALSE)
  if (!all(vapply(mp_list, is.character, logical(1)))) stop("Each gene set must be a character vector of gene symbols.", call. = FALSE)
  mp_list <- lapply(mp_list, function(g) unique(g[!is.na(g) & g != ""]))
  small <- vapply(mp_list, length, integer(1)) < min_genes
  if (any(small)) stop("Gene sets with fewer than ", min_genes, " genes: ", paste(nm[small], collapse = ", "), call. = FALSE)
  if (!is.null(genes)) {
    cov <- vapply(mp_list, function(g) mean(g %in% genes), numeric(1))
    if (any(cov < 0.7)) warning("Low gene coverage for MPs: ", paste(names(cov)[cov < 0.7], round(cov[cov < 0.7], 2), collapse = ", "), call. = FALSE)
  }
  invisible(mp_list)
}

#' Add MP scores to an epithelial Seurat object
#'
#' Runs `Seurat::AddModuleScore` for each MP and stores the scores in metadata
#' columns named after the MPs (e.g. `MP_1_Keratinization`). The object must be
#' log-normalised (`Seurat::NormalizeData`).
#' @param epi Seurat object of epithelial cells.
#' @param mp_list Named list of MP gene sets.
#' @param assay Seurat assay.
#' @export
compute_MP_scores <- function(epi, mp_list = default_MP_list(), assay = "RNA") {
  mp_list <- validate_MP_list(mp_list, genes = rownames(epi))
  tmp <- "sysHIXtmpmp"
  epi <- Seurat::AddModuleScore(epi, features = mp_list, name = tmp, assay = assay)
  old <- paste0(tmp, seq_along(mp_list))
  epi@meta.data[, names(mp_list)] <- epi@meta.data[, old]
  epi@meta.data[, old] <- NULL
  epi
}

#' Check and summarise input before scoring
#'
#' @param epi Epithelial Seurat object (log-normalised).
#' @param tme Non-epithelial Seurat object or metadata with a cell-state column; may be `NULL`
#'   (then no M layer or sysHIX is produced).
#' @param sample_col Metadata column with sample id.
#' @param tissue_col Metadata column with tissue label.
#' @param state_col Metadata column with cell-state label.
#' @param mode Analysis mode.
#' @param nt_label Label of normal tissue in the tissue column.
#' @return Invisibly, a list of per-sample cell counts; messages and warnings flag problems.
#' @export
check_input <- function(epi, tme = NULL, sample_col = "sample_id", tissue_col = "tissue",
                        state_col = "cell_state", mode = "pretrained", nt_label = "NT") {
  md <- .get_meta(epi)
  if (!sample_col %in% names(md)) stop("Epithelial metadata lacks sample column '", sample_col, "'.")
  if (mode != "pretrained") {
    if (!tissue_col %in% names(md)) stop("Mode '", mode, "' needs the tissue column '", tissue_col, "'.")
    if (!nt_label %in% md[[tissue_col]]) stop("No cells labelled '", nt_label, "' in '", tissue_col, "'.")
  }
  res <- list(epi = table(md[[sample_col]]))
  if (!is.null(tme)) {
    tm <- .get_meta(tme)
    if (!all(c(sample_col, state_col) %in% names(tm)))
      stop("Non-epithelial metadata needs columns '", sample_col, "' and '", state_col, "'.")
    res$tme <- table(tm[[sample_col]][!is.na(tm[[state_col]])])
  }
  invisible(res)
}
