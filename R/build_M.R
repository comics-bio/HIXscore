# Mode 3, M layer: build cell-state composition ecotypes (CMs) with NMF

#' Sample x cell-state composition matrix for M-layer training
#'
#' Proportions of each immune/stromal cell state per sample (denominator: cells in the
#' retained states). States are kept if, in at least one tissue group, their mean proportion
#' is >= `min_mean_prop` and they occur in >= `min_prevalence` of that group's samples
#' (the rule used for the shipped model). Without `tissue_col`, all samples form one group.
#' @param x Seurat object or cell-level metadata of non-epithelial cells.
#' @param sample_col Metadata column with sample id.
#' @param state_col Metadata column with cell-state label.
#' @param tissue_col Metadata column with tissue label.
#' @param min_cells Minimum cells per sample.
#' @param min_mean_prop Minimum mean proportion of a state.
#' @param min_prevalence Minimum fraction of samples containing a state.
#' @return States x samples matrix of proportions.
#' @export
build_M_composition <- function(x, sample_col = "sample_id", state_col = "cell_state",
                                tissue_col = NULL, min_cells = 200,
                                min_mean_prop = 0.01, min_prevalence = 0.10) {
  md <- .get_meta(x)
  md <- md[!is.na(md[[state_col]]), , drop = FALSE]
  s <- as.character(md[[sample_col]]); st <- sub("^Fibroblasts_", "Fibro_", as.character(md[[state_col]]))
  tot <- table(s); ok <- s %in% names(tot)[tot >= min_cells]
  if (!any(ok)) stop("No sample has at least ", min_cells, " immune/stromal cells.")
  s <- s[ok]; st <- st[ok]
  grp <- if (is.null(tissue_col)) rep("all", length(s)) else as.character(md[[tissue_col]])[ok]
  P <- sweep(unclass(table(st, s)), 2, colSums(unclass(table(st, s))), "/")
  g_of_s <- grp[match(colnames(P), s)]
  keep <- rep(FALSE, nrow(P))
  for (g in unique(g_of_s)) {
    Pg <- P[, g_of_s == g, drop = FALSE]
    keep <- keep | (rowMeans(Pg) >= min_mean_prop & rowMeans(Pg > 0) >= min_prevalence)
  }
  if (sum(keep) < 5) stop("Fewer than 5 cell states pass the filters; relax min_mean_prop / min_prevalence.")
  P <- P[keep, , drop = FALSE]
  sweep(P, 2, colSums(P), "/")
}

.minmax_params <- function(P) list(min = apply(P, 1, min), max = apply(P, 1, max))
.minmax_scale <- function(P, mm) {
  rng <- pmax(mm$max - mm$min, 1e-12)
  (P - mm$min[rownames(P)]) / rng[rownames(P)]
}

#' Select the optimal NMF rank from cophenetic correlations
#'
#' A rank k is a candidate when the cophenetic correlation rho rises up to k and then falls:
#' rho(k-2) < rho(k-1) < rho(k) and rho(k) > rho(k+1) (strict inequalities). Candidates therefore
#' need two smaller and one larger rank in the scan. If several ranks qualify, the one with the
#' highest rho is chosen (ties: the smaller rank).
#' @param rho Cophenetic correlations, ordered by rank.
#' @param ranks The corresponding ranks.
#' @return list(optimal, candidates); `optimal` is `NA` if no rank qualifies.
#' @export
select_M_rank <- function(rho, ranks) {
  o <- order(ranks); rho <- rho[o]; ranks <- ranks[o]
  n <- length(rho); cand <- logical(n)
  if (n >= 4) for (i in 3:(n - 1))
    cand[i] <- isTRUE(rho[i - 2] < rho[i - 1] && rho[i - 1] < rho[i] && rho[i] > rho[i + 1])
  cands <- ranks[cand]
  opt <- if (length(cands)) cands[which.max(rho[cand])] else NA_integer_
  list(optimal = opt, candidates = cands)
}

#' Scan NMF ranks and select the optimal rank
#'
#' Fits `nsNMF` over `ranks` on the min-max scaled composition and returns the stability/fit
#' metrics. The optimal rank is chosen by [select_M_rank()] from the cophenetic correlation
#' (rho(k-2) < rho(k-1) < rho(k) and rho(k) > rho(k+1)). The scan must therefore cover ranks
#' on both sides of the expected optimum (default 2:10).
#' @param P States x samples proportions from [build_M_composition()].
#' @param ranks Ranks to scan (at least 4).
#' @param nrun Number of NMF runs per rank.
#' @param seed Random seed.
#' @param method NMF algorithm.
#' @return data.frame of metrics per rank with a logical column `candidate`; attributes
#'   `optimal_rank` (NA if none qualifies) and `fit` (the NMF object).
#' @export
scan_M_rank <- function(P, ranks = 2:10, nrun = 50, seed = 123456, method = "nsNMF") {
  .need("NMF", "NMF rank scanning")
  if (length(ranks) < 4) stop("At least 4 ranks are needed to apply the rank-selection rule.")
  V <- .minmax_scale(P, .minmax_params(P))
  fit <- NMF::nmf(V, rank = ranks, method = method, seed = seed, nrun = nrun, .options = "v-p")
  m <- fit$measures
  out <- data.frame(rank = m$rank, cophenetic = m$cophenetic, dispersion = m$dispersion,
                    silhouette = m$silhouette.consensus, rss = m$rss, row.names = NULL)
  sel <- select_M_rank(out$cophenetic, out$rank)
  out$candidate <- out$rank %in% sel$candidates
  if (is.na(sel$optimal))
    warning("No rank satisfies rho(k-2) < rho(k-1) < rho(k) and rho(k) > rho(k+1); widen `ranks`, increase `nrun`, or choose k manually.", call. = FALSE)
  attr(out, "optimal_rank") <- sel$optimal
  attr(out, "fit") <- fit
  out
}

#' Build an M-layer model (CM ecotypes)
#'
#' Either fits NMF with `k` factors on the min-max scaled composition, or accepts a
#' user-supplied basis matrix `W` (states x k; e.g. from your own NMF) and only derives
#' the frozen scaling. The result has the same structure as [default_M_model()] and
#' can be passed to [compute_M_CM_scores()].
#' @param P Proportion matrix from [build_M_composition()].
#' @param k Number of CMs, or `"auto"` to run [scan_M_rank()] and use its optimal rank (ignored if `W` is given).
#' @param W Optional user basis matrix with state names as rownames.
#' @param nrun Number of NMF runs.
#' @param seed Random seed.
#' @param method Algorithm.
#' @export
build_M_ecotypes <- function(P, k = NULL, W = NULL, nrun = 100, seed = 1234, method = "nsNMF") {
  mm <- .minmax_params(P); V <- .minmax_scale(P, mm)
  if (is.null(W)) {
    if (is.null(k)) stop("Give k (see scan_M_rank()), k = \"auto\", or a user-supplied W.")
    if (identical(k, "auto")) {
      sc <- scan_M_rank(P, nrun = max(nrun %/% 2, 20), seed = seed, method = method)
      k <- attr(sc, "optimal_rank")
      if (is.na(k)) stop("Automatic rank selection found no rank; inspect scan_M_rank() and set k manually.")
      message("Selected rank k = ", k, " (cophenetic criterion).")
    }
    .need("NMF", "NMF ecotype building")
    fit <- NMF::nmf(V, k, method = method, seed = seed, nrun = nrun, .options = "-p")
    W <- NMF::basis(fit)
  } else {
    W <- as.matrix(W)
    if (is.null(rownames(W))) stop("W needs state names as rownames.")
    miss <- setdiff(rownames(W), rownames(P))
    if (length(miss)) stop("W contains states absent from the composition matrix: ", paste(utils::head(miss), collapse = ", "))
    if (any(W < 0)) stop("W must be non-negative.")
    P <- P[rownames(W), , drop = FALSE]; P <- sweep(P, 2, colSums(P), "/")
    mm <- .minmax_params(P); V <- .minmax_scale(P, mm)
  }
  rownames(W) <- rownames(V); colnames(W) <- paste0("CM", seq_len(ncol(W)))
  H <- vapply(seq_len(ncol(V)), function(i) { h <- nnls::nnls(W, V[, i])$x; if (sum(h) > 0) h / sum(h) else rep(NA_real_, ncol(W)) }, numeric(ncol(W)))
  H <- matrix(H, nrow = ncol(W), dimnames = list(colnames(W), colnames(V)))
  list(W = W, comp_min = mm$min, comp_max = mm$max, k = ncol(W), V_train = V, H_train = H,
       min_total_cells = NA_integer_, creation_date = Sys.Date())
}
