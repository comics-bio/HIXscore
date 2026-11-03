# Entry points: MP scoring, input checks and the one-call pipeline

#' Epithelial meta-program (MP) gene lists shipped with the package
#' @export
default_MP_list <- function() {
  readRDS(system.file("extdata", "MP_list.rds", package = "HIXscore", mustWork = TRUE))
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
  cov <- vapply(mp_list, function(g) mean(g %in% rownames(epi)), numeric(1))
  if (any(cov < 0.7))
    warning("Low gene coverage for MPs: ", paste(names(cov)[cov < 0.7], round(cov[cov < 0.7], 2), collapse = ", "), call. = FALSE)
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
