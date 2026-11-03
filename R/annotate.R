# Cell annotation by two-level Seurat label transfer (port of scRNA_cell_Transfer_*.R).
# Level 1: major cell type per cell; level 2: cell state within each major type.
# Query samples are annotated one at a time, so a sample's labels never depend on other samples.

.drop_noise_genes <- function(x) {
  cc <- Seurat::cc.genes.updated.2019
  x[!x %in% c(cc$s.genes, cc$g2m.genes) & !grepl("^MT-|^MT\\.|^LOC(0-9)", x)]
}

.prep_norm <- function(obj, nfeat = 3000) {
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  Seurat::FindVariableFeatures(obj, nfeatures = min(nfeat, nrow(obj)), verbose = FALSE)
}

# TransferData with a k.weight back-off (k.weight is reduced by a factor 0.6 until the transfer succeeds)
.transfer_retry <- function(anchors, refdata, dims, k.weight, min_k) {
  pred <- NULL
  while (k.weight >= min_k && is.null(pred)) {
    pred <- tryCatch(Seurat::TransferData(anchorset = anchors, refdata = refdata, dims = dims,
                                          k.weight = k.weight, verbose = FALSE),
                     error = function(e) NULL)
    if (is.null(pred)) k.weight <- floor(k.weight * 0.6)
  }
  pred
}

#' Prepare a full annotation reference
#'
#' Splits the reference by major cell type and normalises each subset once, so that annotating many
#' samples does not repeat this work.
#'
#' A reference of your own has to contain:
#' * a Seurat (v4) object with the raw counts in the `RNA` assay (all cells, all genes) and gene symbols as
#'   row names;
#' * metadata columns with the major cell type (`major_col`, default `celltype_major`) and the cell state
#'   (`minor_col`, default `cell_type2`), without missing values; epithelial cells must be labelled
#'   `"Epithelial"` in the major-type column, and the cell states must be nested in the major types;
#' * a PCA reduction (`reduction`, default `"pca"`) with the embeddings of **all** cells and the **feature
#'   loadings** (the genes used for the PCA, usually the variable genes);
#' * at least `min_cells_per_type` (default 50) cells per cell state, and at least 100 cells per major type;
#'   smaller states and types are not used for the second level.
#'
#' For the shipped M-layer model (modes 1 and 2) the cell states must include the 72 immune and stromal states
#' of `rownames(default_M_model()$W)` ([check_reference()] lists the missing ones). A reference with a different
#' set of cell states can only be used when you build your own M layer (see [hixscore_train()]).
#'
#' With `reduction = "pca"` (default) the second level projects query cells onto a PCA that is
#' recomputed within each major type: each subset is additionally scaled and reduced by PCA with
#' `min(50, half the cells, genes - 1)` (at least 10) components.
#'
#' @param ref Seurat reference object.
#' @param ref_features Named list (names = major types) of marker genes added to the transfer
#'   features of the second level (e.g. top markers of each cell state).
#' @param major_col,minor_col Metadata columns with major type and cell state.
#' @param reduction Reduction of `ref` used to project query cells (default `"pca"`).
#' @param min_cells_per_type Minimum reference cells of a cell state.
#' @return Object of class `HIXscore_reference`.
#' @export
prepare_reference <- function(ref, ref_features = NULL, major_col = "celltype_major",
                              minor_col = "cell_type2", reduction = "pca", min_cells_per_type = 50) {
  stopifnot(all(c(major_col, minor_col) %in% colnames(ref@meta.data)), reduction %in% names(ref@reductions))
  Seurat::DefaultAssay(ref) <- "RNA"
  parts <- Seurat::SplitObject(ref, split.by = major_col)
  minor <- lapply(parts, function(p) {
    cnt <- table(p@meta.data[[minor_col]]); ok <- names(cnt)[cnt >= min_cells_per_type]
    if (length(ok) == 0 || ncol(p) < 100) return(NULL)
    p <- p[, p@meta.data[[minor_col]] %in% ok]
    if (length(ok) == 1) return(list(single = ok))
    p <- .prep_norm(p, 2000)
    if (identical(reduction, "pca")) p <- .prep_pca(p)
    p
  })
  structure(list(major = ref, minor = minor, ref_features = ref_features, major_col = major_col,
                 minor_col = minor_col, reduction = reduction,
                 states = sort(unique(ref@meta.data[[minor_col]])), info = list(created = Sys.time())),
            class = "HIXscore_reference")
}

# scale + PCA within a major-type subset (level-2 reference preprocessing)
.prep_pca <- function(p) {
  npcs <- max(min(50, floor(ncol(p) * 0.5), nrow(p) - 1), 10)
  p <- Seurat::ScaleData(p, verbose = FALSE)
  Seurat::RunPCA(p, npcs = npcs, verbose = FALSE)
}

.transfer_major <- function(ref, que, major_col, reduction, dims = 1:30, k.weight = 50) {
  feat <- unique(c(Seurat::VariableFeatures(ref), Seurat::VariableFeatures(que)))
  feat <- .drop_noise_genes(intersect(intersect(rownames(ref), rownames(que)), feat))
  anchors <- Seurat::FindTransferAnchors(reference = ref, query = que, features = feat,
               normalization.method = "LogNormalize", reference.reduction = reduction, dims = dims, verbose = FALSE)
  pred <- .transfer_retry(anchors, ref@meta.data[[major_col]], dims, k.weight, 10)
  if (is.null(pred)) return(NULL)
  data.frame(predicted.major = pred$predicted.id, prediction.score.major = pred$prediction.score.max,
             row.names = rownames(pred), stringsAsFactors = FALSE)
}

# level 2 inside one major type; returns data.frame(predicted.minor, score) or NULL when skipped
.transfer_minor_one <- function(rsub, que_sub, ref_markers, minor_col, reduction, dims = 1:30,
                                min_query_cells = 50, k.filter = NULL, k.weight = NULL) {
  n_que <- ncol(que_sub)
  if (is.null(rsub) || n_que < min_query_cells) return(NULL)
  if (!inherits(rsub, "Seurat"))      # major type with a single cell state: no transfer needed
    return(data.frame(predicted.minor = rsub$single, prediction.score.minor = 1, row.names = colnames(que_sub)))
  que_sub <- .prep_norm(que_sub, 2000)
  max_npcs <- function(o) max(min(50, floor(ncol(o) * 0.5), nrow(o) - 1), 10)
  dims_use <- 1:min(length(dims), max_npcs(rsub), max_npcs(que_sub))
  kf <- if (is.null(k.filter)) min(200, floor(n_que * 0.8)) else min(k.filter, floor(n_que * 0.8))
  common <- intersect(rownames(rsub), rownames(que_sub))
  feat <- intersect(Seurat::VariableFeatures(rsub), Seurat::VariableFeatures(que_sub))
  feat <- .drop_noise_genes(intersect(common, union(feat, ref_markers)))
  anchors <- tryCatch(Seurat::FindTransferAnchors(reference = rsub, features = feat, query = que_sub,
               normalization.method = "LogNormalize", reference.reduction = reduction, dims = dims_use,
               k.filter = kf, verbose = FALSE), error = function(e) NULL)
  if (is.null(anchors)) return(NULL)
  n_anchors <- nrow(anchors@anchors)
  if (n_anchors < 10) return(NULL)
  kw <- if (is.null(k.weight)) min(50, floor(n_anchors * 0.5)) else min(k.weight, floor(n_anchors * 0.5))
  pred <- .transfer_retry(anchors, rsub@meta.data[[minor_col]], dims_use, kw, 5)
  if (is.null(pred)) return(NULL)
  pred <- pred[colnames(que_sub), , drop = FALSE]
  data.frame(predicted.minor = pred$predicted.id, prediction.score.minor = pred$prediction.score.max,
             row.names = rownames(pred), stringsAsFactors = FALSE)
}

#' Annotate cells by label transfer from a reference
#'
#' Two-level annotation as used for the HIXscore reference: (1) major cell type, (2) cell state within
#' each major type. Samples are annotated independently. Cells of the second level that cannot be
#' transferred (too few cells of that major type in the sample, too few anchors) get `NA`.
#'
#' @param obj Seurat object with raw counts.
#' @param reference A `HIXscore_reference` ([prepare_reference()] or [load_reference()]).
#' @param sample_col Metadata column with the sample id (annotation is run per sample).
#' @param min_cells_per_sample Samples with fewer cells are not annotated.
#' @param minor_for Major types for which the cell state is predicted: `"all"` or a character vector.
#'   The default `"non-epithelial"` skips the epithelial states, which no layer of the score uses.
#' @param min_query_cells Minimum cells of a major type in a sample for the second level.
#' @param dims Dimensions of the reference embedding.
#' @param reduction Reduction of the reference used to project query cells; `NULL` (default) uses the
#'   reduction the reference was prepared with ([prepare_reference()], default `"pca"`). A different
#'   reduction requires a reference prepared with it.
#' @param min_score_major,min_score_minor Optional thresholds; cells below are labelled `Low_Quality`
#'   (major) or `<major>_Unassigned` (state). Default `NULL`: no filtering.
#' @param checkpoint_dir Optional directory. The result of each sample is written there as soon as it is finished,
#'   and samples already present are read from it instead of being annotated again, so an interrupted or partly
#'   failed run can be repeated at the cost of the missing samples only. Use a new directory for a different
#'   reference or settings.
#' @param verbose Print progress.
#' @return `obj` with metadata `predicted.major`, `prediction.score.major`, `predicted.minor`,
#'   `prediction.score.minor`, `annotation.confidence`; attributes `annotation_summary` and `annotation_failed`
#'   (samples that could not be annotated and the reason). A sample that fails does not stop the others; its cells
#'   keep `NA` labels.
#' @export
annotate_cells <- function(obj, reference, sample_col = "sample_id", min_cells_per_sample = 50,
                           minor_for = "non-epithelial", min_query_cells = 50, dims = 1:30, reduction = NULL,
                           min_score_major = NULL, min_score_minor = NULL, checkpoint_dir = NULL, verbose = TRUE) {
  .check_seurat()
  stopifnot(inherits(reference, "HIXscore_reference"), sample_col %in% colnames(obj@meta.data))
  if (!is.null(reduction)) {
    if (!reduction %in% names(reference$major@reductions))
      stop("Reference has no reduction '", reduction, "' (available: ", paste(names(reference$major@reductions), collapse = ", "), ").", call. = FALSE)
    if (identical(reduction, "pca") && !identical(reference$reduction, "pca"))
      stop("This reference was prepared with reduction = '", reference$reduction, "', so its level-2 subsets have no ",
           "within-type PCA. Re-run prepare_reference(..., reduction = \"pca\") to annotate with 'pca'.", call. = FALSE)
    reference$reduction <- reduction
  }
  Seurat::DefaultAssay(obj) <- "RNA"
  samp <- as.character(obj@meta.data[[sample_col]]); ids <- unique(samp)
  major_all <- names(reference$minor)
  do_minor <- if (identical(minor_for, "all")) major_all
              else if (identical(minor_for, "non-epithelial")) setdiff(major_all, "Epithelial") else minor_for
  res <- vector("list", length(ids)); names(res) <- ids; failed <- list()
  if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)
  ckpt <- function(s) file.path(checkpoint_dir, paste0(make.names(s), "_", substr(as.character(utils::packageVersion("HIXscore")), 1, 12), ".rds"))
  annotate_one <- function(s, cells) {
    que <- .prep_norm(obj[, cells], 3000)
    mj <- .transfer_major(reference$major, que, reference$major_col, reference$reduction, dims)
    if (is.null(mj)) stop("major-type transfer failed")
    out <- mj[colnames(que), , drop = FALSE]
    out$predicted.minor <- NA_character_; out$prediction.score.minor <- NA_real_
    que$predicted.major <- out$predicted.major
    for (mt in intersect(do_minor, unique(out$predicted.major))) {
      qs <- que[, out$predicted.major == mt]
      mn <- .transfer_minor_one(reference$minor[[mt]], qs, reference$ref_features[[mt]], reference$minor_col,
                                reference$reduction, dims, min_query_cells)
      if (!is.null(mn)) { out[rownames(mn), "predicted.minor"] <- mn$predicted.minor
                          out[rownames(mn), "prediction.score.minor"] <- mn$prediction.score.minor }
    }
    out
  }
  for (s in ids) {
    cells <- colnames(obj)[samp == s]
    if (length(cells) < min_cells_per_sample) {
      failed[[s]] <- paste0(length(cells), " cells < min_cells_per_sample (", min_cells_per_sample, ")")
      if (verbose) message("[", s, "] ", length(cells), " cells < ", min_cells_per_sample, ": not annotated")
      next
    }
    if (!is.null(checkpoint_dir) && file.exists(ckpt(s))) {
      prev <- readRDS(ckpt(s))
      if (identical(sort(rownames(prev)), sort(cells))) { res[[s]] <- prev; if (verbose) message("[", s, "] read from checkpoint"); next }
    }
    out <- tryCatch(annotate_one(s, cells), error = function(e) { failed[[s]] <<- conditionMessage(e); NULL })
    if (is.null(out)) { warning("[", s, "] not annotated: ", failed[[s]], call. = FALSE); next }
    res[[s]] <- out
    if (!is.null(checkpoint_dir)) saveRDS(out, ckpt(s))
    if (verbose) message("[", s, "] ", length(cells), " cells annotated")
  }
  done <- res[!vapply(res, is.null, logical(1))]
  if (length(done) == 0) stop("No sample could be annotated.")
  ann <- do.call(rbind, unname(done))     # row names are the (unique) cell names
  if (anyDuplicated(colnames(obj))) stop("Cell names must be unique across samples.")
  md <- ann[match(colnames(obj), rownames(ann)), , drop = FALSE]; rownames(md) <- colnames(obj)
  md$annotation.confidence <- ifelse(is.na(md$prediction.score.minor), "unassigned",
                               ifelse(md$prediction.score.minor > 0.6, "high", ifelse(md$prediction.score.minor > 0.4, "medium", "low")))
  if (!is.null(min_score_major)) md$predicted.major[!is.na(md$prediction.score.major) & md$prediction.score.major < min_score_major] <- "Low_Quality"
  if (!is.null(min_score_minor)) { lo <- !is.na(md$prediction.score.minor) & md$prediction.score.minor < min_score_minor
    md$predicted.minor[lo] <- paste0(md$predicted.major[lo], "_Unassigned") }
  for (v in colnames(md)) obj[[v]] <- md[[v]]
  attr(obj, "annotation_summary") <- annotation_summary(obj)
  attr(obj, "annotation_failed") <- if (length(failed)) data.frame(sample = names(failed), reason = unlist(failed), row.names = NULL) else NULL
  obj
}

#' Summarise an annotation
#' @param obj Annotated Seurat object or metadata.
#' @param model M-layer model whose cell states define the vocabulary (default shipped model).
#' @return list(n_cells, annotated_fraction, major_counts, in_model_state_fraction).
#' @export
annotation_summary <- function(obj, model = default_M_model()) {
  md <- .get_meta(obj); a <- !is.na(md$predicted.major)
  nonepi <- a & md$predicted.major != "Epithelial"
  list(n_cells = nrow(md), annotated_fraction = mean(a), major_counts = table(md$predicted.major),
       in_model_state_fraction = if (any(nonepi)) mean(md$predicted.minor[nonepi] %in% rownames(model$W), na.rm = FALSE) else NA_real_)
}
