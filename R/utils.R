# Internal helpers: supported Seurat version and optional-package checks

# HIXscore supports Seurat 4.x (4.3 - 4.4). Seurat >= 5 changed assay classes and the
# slot/layer interface, so we stop early instead of returning possibly wrong scores.
.check_seurat <- function() {
  v <- tryCatch(utils::packageVersion("SeuratObject"), error = function(e) NULL)
  s <- tryCatch(utils::packageVersion("Seurat"), error = function(e) NULL)
  if (is.null(s) || is.null(v)) stop("Package 'Seurat' (4.x) is required.", call. = FALSE)
  if (s >= "5.0.0" || v >= "5.0.0")
    stop("HIXscore supports Seurat 4.x only (4.3 - 4.4); found Seurat ", s,
         " / SeuratObject ", v, ". Please use a Seurat 4 environment.", call. = FALSE)
  invisible(TRUE)
}

# error early, with an actionable message, when an optional package is missing
.need <- function(pkgs, why = NULL) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop("Package", if (length(miss) > 1) "s " else " ", paste0("'", miss, "'", collapse = ", "),
         " required", if (!is.null(why)) paste0(" for ", why), ". Install with install.packages() or BiocManager::install().",
         call. = FALSE)
  invisible(TRUE)
}
