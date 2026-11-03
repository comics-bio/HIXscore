# Mode 3, C layer: derive baseline-deviation module (BDM) gene programs.
# Re-implements step_3_baseline_MPs/01HAPs_select.R as parameterised functions:
#   1. pseudobulk per sample x epithelial state
#   2. leave-one-cohort-out (LOCO) variance partitioning on NT -> "constrained" genes
#   3. LOCO disease-vs-NT deviation -> consistently up/down genes
#   4. co-deviation clustering (consensus clustering) -> modules

#' Pseudobulk expression per sample x epithelial state
#'
#' @param epi Seurat object of epithelial cells (raw counts available).
#' @param state_col Epithelial state label per cell (e.g. the MP state); unresolved cells should be removed beforehand.
#' @param min_cells Minimum cells per sample x state stratum (default 10).
#' @param min_logcpm Genes with mean log2(CPM+1) below this are dropped.
#' @param exclude_patterns Regex patterns of technical genes to drop.
#' @param sample_col Metadata column with sample id.
#' @param tissue_col Metadata column with tissue label.
#' @param series_col Metadata column with cohort/study.
#' @param subsite_col Optional metadata column with anatomical subsite.
#' @param assay Seurat assay.
#' @return list(counts, meta, logcpm) with columns named `<sample>_<state>`.
#' @export
bdm_pseudobulk <- function(epi, sample_col = "sample_id", state_col = "MP_state", tissue_col = "tissue",
                           series_col = "series", subsite_col = NULL, min_cells = 10, min_logcpm = 1,
                           exclude_patterns = c("^RPL", "^RPS", "^MT-", "^HBA", "^HBB", "^MALAT1", "^NEAT1", "^XIST", "^ENSG", "^LINC"),
                           assay = "RNA") {
  .need(c("edgeR", "limma", "Matrix"))
  md <- epi@meta.data
  need <- c(sample_col, state_col, tissue_col, series_col, subsite_col)
  if (!all(need %in% names(md))) stop("Missing metadata columns: ", paste(setdiff(need, names(md)), collapse = ", "))
  key <- paste(md[[sample_col]], md[[state_col]], sep = "_")
  n <- table(key); ok <- key %in% names(n)[n >= min_cells] & !is.na(md[[state_col]])
  cnt <- suppressWarnings(SeuratObject::GetAssayData(epi, assay = assay, slot = "counts"))[, ok, drop = FALSE]
  key <- key[ok]; md <- md[ok, , drop = FALSE]
  f <- factor(key); G <- Matrix::sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), nlevels(f)),
                                              dimnames = list(NULL, levels(f)))
  pb <- as.matrix(cnt %*% G)
  first <- !duplicated(key)
  meta <- data.frame(pb_id = key[first], sample_id = as.character(md[[sample_col]][first]),
                     MP_state = as.character(md[[state_col]][first]), tissue = as.character(md[[tissue_col]][first]),
                     series = as.character(md[[series_col]][first]),
                     subsite = if (is.null(subsite_col)) "all" else as.character(md[[subsite_col]][first]),
                     n_cells = as.integer(table(key)[key[first]]), stringsAsFactors = FALSE)
  pb <- pb[, meta$pb_id, drop = FALSE]
  tech <- grepl(paste(exclude_patterns, collapse = "|"), rownames(pb), ignore.case = TRUE)
  low <- rowMeans(log2(edgeR::cpm(pb) + 1)) < min_logcpm
  pb <- pb[!tech & !low, , drop = FALSE]
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = pb), method = "TMM")
  list(counts = pb, meta = meta, logcpm = limma::voom(dge, plot = FALSE)$E)
}

# one LOCO iteration set: pseudobulk is recomputed on the training cohorts minus the left-out one
.loco_sets <- function(meta, min_samples = 5, loco_cohorts = NULL) {
  per <- tapply(meta$sample_id, meta$series, function(x) length(unique(x)))
  cand <- if (is.null(loco_cohorts)) names(per)[per > min_samples] else loco_cohorts
  if (length(cand) < 2)
    stop("Leave-one-cohort-out selection needs at least 2 cohorts (`series`) with more than ", min_samples,
         " samples each; found ", length(cand), ". Cohort sizes (samples): ", paste(names(per), per, sep = "=", collapse = ", "), ".", call. = FALSE)
  if (length(cand) < 3) warning("Only ", length(cand), " cohorts available for leave-one-cohort-out selection; the gene frequency across iterations is coarse.", call. = FALSE)
  cand
}

#' Constrained (stable-in-normal-tissue) genes by leave-one-cohort-out variance partitioning
#'
#' For each LOCO iteration, variancePartition is run on NT pseudobulk with random effects for
#' state, sample, subsite and cohort; genes whose variance fraction is below `max_var` for all
#' four factors are "constrained". Genes constrained in >= `min_freq` of the iterations are kept.
#' @param epi Seurat epithelial object (as for [bdm_pseudobulk()]); pseudobulk is rebuilt per iteration.
#' @param train_cohorts Cohorts (`series`) used for training; `NULL` = all.
#' @param nt_label Label of normal tissue.
#' @param max_var Maximum variance fraction explained by state, sample, subsite or cohort.
#' @param min_freq Minimum fraction of LOCO iterations in which a gene must be selected.
#' @param min_nt_pseudobulk Minimum NT pseudobulk samples for an iteration.
#' @param workers Number of parallel workers.
#' @param ... Passed to [bdm_pseudobulk()] (column names, thresholds).
#' @export
bdm_constrained_genes <- function(epi, train_cohorts = NULL, nt_label = "NT", max_var = 0.15, min_freq = 0.7,
                                  min_nt_pseudobulk = 10, workers = 1, ...) {
  .need(c("variancePartition", "BiocParallel"))
  md <- epi@meta.data; series_col <- list(...)$series_col; if (is.null(series_col)) series_col <- "series"
  tr <- if (is.null(train_cohorts)) unique(md[[series_col]]) else train_cohorts
  pb_all <- bdm_pseudobulk(epi[, md[[series_col]] %in% tr], ...)
  sel <- .loco_sets(pb_all$meta)
  bp <- if (workers > 1) BiocParallel::MulticoreParam(workers) else BiocParallel::SerialParam()
  res <- list()
  for (lo in sel) {
    pb <- bdm_pseudobulk(epi[, md[[series_col]] %in% setdiff(tr, lo)], ...)
    m <- pb$meta[pb$meta$tissue == nt_label, ]
    if (nrow(m) < min_nt_pseudobulk) { message("LOCO ", lo, ": too few NT pseudobulk samples, skipped"); next }
    for (v in c("MP_state", "sample_id", "subsite", "series")) m[[v]] <- factor(m[[v]])
    rownames(m) <- m$pb_id
    facs <- c("MP_state", "sample_id", "subsite", "series")
    facs <- facs[vapply(facs, function(v) nlevels(droplevels(m[[v]])) > 1, logical(1))]   # drop single-level factors
    fml <- stats::as.formula(paste("~", paste0("(1 | ", facs, ")", collapse = " + ")))
    vp <- tryCatch(as.data.frame(variancePartition::fitExtractVarPartModel(
      pb$logcpm[, m$pb_id], fml, m, BPPARAM = bp, showWarnings = FALSE)),
      error = function(e) { message("LOCO ", lo, ": variancePartition failed: ", conditionMessage(e)); NULL })
    if (is.null(vp)) next
    ok <- rep(TRUE, nrow(vp)); for (v in facs) ok <- ok & vp[[v]] < max_var
    res[[lo]] <- rownames(vp)[ok]
  }
  if (!length(res)) stop("No LOCO iteration produced constrained genes.")
  all <- unique(unlist(res))
  freq <- vapply(all, function(g) mean(vapply(res, function(x) g %in% x, logical(1))), numeric(1))
  info <- data.frame(gene = all, frequency = freq, selected = freq >= min_freq, row.names = NULL)
  list(genes = info$gene[info$selected], frequency = info[order(-info$frequency), ], by_iteration = res)
}

# disease-vs-NT deviation per pseudobulk (relative to the NT mean of the same state)
.deviation_matrix <- function(logcpm, meta, genes, nt_label) {
  genes <- intersect(genes, rownames(logcpm))
  nt <- meta[meta$tissue == nt_label, ]
  ref <- lapply(split(nt$pb_id, nt$MP_state), function(ids) rowMeans(logcpm[genes, ids, drop = FALSE]))
  dis <- meta[meta$tissue != nt_label & meta$MP_state %in% names(ref), ]
  out <- lapply(split(dis, dis$tissue), function(d)
    vapply(seq_len(nrow(d)), function(i) logcpm[genes, d$pb_id[i]] - ref[[d$MP_state[i]]], numeric(length(genes)),
           USE.NAMES = FALSE) |> `dimnames<-`(list(genes, d$pb_id)))
  out
}

#' Consistently deviating genes by leave-one-cohort-out disease-vs-NT analysis
#'
#' In every LOCO iteration and disease type, a gene is "up" ("down") if its deviation from the matched
#' NT state mean is > log2(`min_fc`) (< -log2) with the same sign in >= `min_consistency` of the pseudobulk
#' samples (and >= 3 samples); it must satisfy this in >= `min_disease_types` disease types, and in >= `min_freq` of the iterations.
#' @param constrained Genes from [bdm_constrained_genes()].
#' @param epi Seurat object of epithelial cells.
#' @param train_cohorts Cohorts (`series`) used for training; `NULL` = all.
#' @param nt_label Label of normal tissue in the tissue column.
#' @param min_fc Minimum fold change (linear scale).
#' @param min_consistency Minimum fraction of pseudobulk samples with the same deviation sign.
#' @param min_disease_types Minimum number of disease types in which a gene must deviate.
#' @param min_freq Minimum fraction of LOCO iterations.
#' @param ... Passed to [bdm_pseudobulk()].
#' @export
bdm_deviation_genes <- function(epi, constrained, train_cohorts = NULL, nt_label = "NT", min_fc = 1.2,
                                min_consistency = 0.7, min_disease_types = 3, min_freq = 0.7, ...) {
  md <- epi@meta.data; series_col <- list(...)$series_col; if (is.null(series_col)) series_col <- "series"
  tr <- if (is.null(train_cohorts)) unique(md[[series_col]]) else train_cohorts
  sel <- .loco_sets(bdm_pseudobulk(epi[, md[[series_col]] %in% tr], ...)$meta)
  thr <- log2(min_fc); iters <- list()
  for (lo in sel) {
    pb <- bdm_pseudobulk(epi[, md[[series_col]] %in% setdiff(tr, lo)], ...)
    dev <- .deviation_matrix(pb$logcpm, pb$meta, constrained, nt_label)
    up <- lapply(dev, function(D) { ok <- rowSums(!is.na(D)) >= 3
      rownames(D)[ok & rowMeans(D > 0, na.rm = TRUE) >= min_consistency & rowMeans(D, na.rm = TRUE) > thr] })
    dn <- lapply(dev, function(D) { ok <- rowSums(!is.na(D)) >= 3
      rownames(D)[ok & rowMeans(D < 0, na.rm = TRUE) >= min_consistency & rowMeans(D, na.rm = TRUE) < -thr] })
    pick <- function(l) { t <- table(unlist(l)); names(t)[t >= min(min_disease_types, length(l))] }
    iters[[lo]] <- list(up = pick(up), down = pick(dn))
  }
  fr <- function(k) { g <- unique(unlist(lapply(iters, `[[`, k)))
    f <- vapply(g, function(x) mean(vapply(iters, function(i) x %in% i[[k]], logical(1))), numeric(1)); f[f >= min_freq] }
  up <- fr("up"); dn <- fr("down")
  list(genes = c(names(up), names(dn)), direction = c(stats::setNames(rep("UP", length(up)), names(up)),
                                                      stats::setNames(rep("DOWN", length(dn)), names(dn))),
       frequency = c(up, dn))
}

#' Scan consensus-clustering k for the co-deviation modules
#' @param dev Genes x samples deviation matrix (all disease types), e.g. `bdm_full_deviation()`.
#' @param maxK Maximum number of clusters.
#' @param reps Number of resamplings.
#' @param seed Random seed.
#' @param plot Write consensus plots.
#' @export
bdm_scan_modules <- function(dev, maxK = 10, reps = 200, seed = 42, plot = FALSE) {
  .need("ConsensusClusterPlus")
  cm <- stats::cor(t(dev), use = "pairwise.complete.obs"); cm[is.na(cm)] <- 0
  cc <- ConsensusClusterPlus::ConsensusClusterPlus(as.matrix(stats::as.dist(1 - cm)), maxK = maxK, reps = reps,
         pItem = 0.8, pFeature = 1, clusterAlg = "hc", distance = "euclidean", seed = seed, plot = if (plot) "png" else NULL)
  # mean consensus per cluster as a stability summary
  stab <- vapply(2:maxK, function(k) {
    cmk <- cc[[k]]$consensusMatrix; cl <- cc[[k]]$consensusClass
    mean(vapply(unique(cl), function(c) mean(cmk[cl == c, cl == c]), numeric(1))) }, numeric(1))
  list(consensus = cc, stability = data.frame(k = 2:maxK, mean_within_consensus = stab))
}

#' Deviation matrix on the full training set (all disease types pooled)
#' @param epi Seurat object of epithelial cells.
#' @param genes Genes to include.
#' @param train_cohorts Cohorts (`series`) used for training; `NULL` = all.
#' @param nt_label Label of normal tissue in the tissue column.
#' @param max_na Maximum fraction of missing deviations per gene.
#' @param ... Passed to [bdm_pseudobulk()].
#' @export
bdm_full_deviation <- function(epi, genes, train_cohorts = NULL, nt_label = "NT", max_na = 0.5, ...) {
  md <- epi@meta.data; series_col <- list(...)$series_col; if (is.null(series_col)) series_col <- "series"
  tr <- if (is.null(train_cohorts)) unique(md[[series_col]]) else train_cohorts
  pb <- bdm_pseudobulk(epi[, md[[series_col]] %in% tr], ...)
  D <- do.call(cbind, .deviation_matrix(pb$logcpm, pb$meta, genes, nt_label))
  D[rowMeans(is.na(D)) < max_na, , drop = FALSE]
}

#' Build the BDM gene programs
#'
#' Clusters genes by co-deviation (consensus clustering, `k` modules), assigns each module the
#' majority direction of its genes and keeps only genes with that direction.
#' @param dev Deviation matrix ([bdm_full_deviation()]).
#' @param direction Named UP/DOWN vector from [bdm_deviation_genes()].
#' @param k Number of modules (see [bdm_scan_modules()]).
#' @param drop_modules Cluster ids to discard (e.g. noisy clusters), after inspection.
#' @param reps Number of resamplings.
#' @param seed Random seed.
#' @param scan Optional result of [bdm_scan_modules()] on the same `dev`; its consensus clustering is reused
#'   instead of being recomputed (`k` must not exceed the `maxK` of the scan).
#' @return Named list `Module_1..` of gene vectors, usable as `programs` in [fit_C_baseline()].
#' @export
build_BDM_geneset <- function(dev, direction, k, drop_modules = integer(), reps = 1000, seed = 42, scan = NULL) {
  .need("ConsensusClusterPlus")
  if (is.null(scan)) scan <- bdm_scan_modules(dev, maxK = max(k, 3), reps = reps, seed = seed)
  if (k > length(scan$consensus) - 1) stop("`k` = ", k, " is larger than the maxK of the scan (", length(scan$consensus) - 1, ").", call. = FALSE)
  cl <- scan$consensus[[k]]$consensusClass
  dr <- direction[names(cl)]
  keep <- list()
  for (m in setdiff(sort(unique(cl)), drop_modules)) {
    g <- names(cl)[cl == m]; main <- names(which.max(table(dr[g])))
    keep[[as.character(m)]] <- g[dr[g] == main]
  }
  names(keep) <- paste0("Module_", seq_along(keep))
  attr(keep, "cluster_ids") <- setdiff(sort(unique(cl)), drop_modules)
  keep
}
