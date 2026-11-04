# hixscore_train(): train the pieces of a user model step by step or all at once.
# Steps (each uses only its own inputs, results are kept and can be passed back as `previous`):
#   M_scan    NMF rank scan of the cell-state composition                (derives the M layer)
#   M         build the CM ecotypes (NMF basis) at rank M_k
#   BDM_scan  leave-one-cohort-out constrained/deviating genes + consensus-clustering scan
#   BDM       gene programs from the scan at k = bdm_k
#   E         E-layer baseline (NT samples)
#   C         C-layer baseline and scHIX thresholds (NT cells)
#   fusion    layer scores, out-of-fold S_C/S_M and the fusion model

.TRAIN_STEPS <- c("M_scan", "M", "BDM_scan", "BDM", "E", "C", "fusion")

.norm_steps <- function(steps) {
  steps <- unique(as.character(steps))
  if ("ALL" %in% toupper(steps)) return(.TRAIN_STEPS)
  bad <- setdiff(steps, .TRAIN_STEPS)
  if (length(bad)) stop("Unknown step(s): ", paste(bad, collapse = ", "), ". Use any of ", paste(.TRAIN_STEPS, collapse = ", "), ", or \"ALL\".", call. = FALSE)
  .TRAIN_STEPS[.TRAIN_STEPS %in% steps]
}

#' Train the layers and the fusion model on your own data
#'
#' Runs the training steps of the model one at a time or together. Steps: `"M_scan"` (NMF rank scan),
#' `"M"` (cell-state ecotypes at rank `M_k`), `"BDM_scan"` (leave-one-cohort-out constrained and deviating
#' genes and the consensus-clustering scan), `"BDM"` (gene programs at `bdm_k`), `"E"` (E-layer baseline),
#' `"C"` (C-layer baseline and scHIX thresholds) and `"fusion"` (layer scores, out-of-fold S_C and S_M, and the
#' fusion model). `"ALL"` runs every step that is needed.
#'
#' `M_model` and `programs` decide whether the M-layer ecotypes and the BDM gene programs are derived from your
#' data (`"derive"`), taken from the package (`"default"`) or supplied by you (an object from [build_M_ecotypes()]
#' or a named list of gene sets). Only `"derive"` runs the corresponding steps. With `"derive"`, the BDM number of
#' modules `bdm_k` has to be chosen from the `BDM_scan` stability table and must be given before the `BDM` step.
#'
#' A step that fails does not stop the others (`on_error = "continue"`); the error is recorded in `status` and the
#' steps that depend on it are skipped. Pass the returned object as `previous` to run further steps or to repeat a
#' failed one without recomputing the finished steps: with `previous`, `"ALL"` runs only the steps that are not
#' finished yet, while steps named explicitly are always run again.
#'
#' @param epi Epithelial Seurat object (log-normalised) with the sample, tissue and series columns. For the BDM
#'   steps it also needs raw counts and the epithelial state label `bdm_state_col`.
#' @param tme Non-epithelial cells (Seurat object or data.frame) with the sample column and `state_col`.
#' @param labels data.frame with `sample_id` and `label` (classes in `levels`), for the `"fusion"` step.
#' @param steps Steps to run: any of `"M_scan"`, `"M"`, `"BDM_scan"`, `"BDM"`, `"E"`, `"C"`, `"fusion"`, or `"ALL"`.
#' @param levels Ordered classes, first = NT. Length 3 (e.g. NT, non-malignant, malignant) or 2 (NT, disease).
#' @param M_model `"default"`, `"derive"`, or an M-layer model object.
#' @param programs `"default"`, `"derive"`, or a named list of gene sets.
#' @param M_k Number of CMs for the `"M"` step: an integer or `"auto"` (cophenetic rule of [select_M_rank()]).
#' @param bdm_k Number of consensus-clustering modules for the `"BDM"` step.
#' @param bdm_drop Cluster ids to discard in the `"BDM"` step (see [build_BDM_geneset()]).
#' @param bdm_state_col Metadata column of `epi` with the epithelial state of each cell (used by the BDM steps).
#' @param bdm_args Named list of further arguments for [bdm_constrained_genes()], [bdm_deviation_genes()] and
#'   [bdm_full_deviation()] (for example `list(min_disease_types = 2)`); `bdm_maxK` and `bdm_reps` set the scan.
#' @param bdm_maxK,bdm_reps Maximum number of clusters and resamplings of the consensus-clustering scan.
#' @param M_ranks Ranks scanned in the `"M_scan"` step.
#' @param previous Result of an earlier `hixscore_train()` call whose step results are reused.
#' @param on_error `"continue"` (default) or `"stop"`.
#' @param sample_col,tissue_col,series_col,state_col Metadata column names.
#' @param nt_label Label of normal tissue.
#' @param min_cells_epi,min_cells_tme Minimum epithelial / immune+stromal cells per sample.
#' @param add_mp_scores Compute MP scores on `epi`.
#' @param M_W Optional user NMF basis matrix of the M layer (cell states x CMs, row names = cell states, non-negative;
#'   see [build_M_ecotypes()]). When given, the `"M"` step builds the M-layer model from this `W` (no rank scan, `M_k`
#'   is ignored) and `"M_scan"` is not run; `tme` is still needed to build the composition matrix.
#' @param mp_list Optional user E-layer gene sets (named list, see [validate_MP_list()]) used instead of
#'   [default_MP_list()] in the `"E"` and `"fusion"` steps; the E baseline is fitted on these sets and the sets are
#'   returned as `mp_list`.
#' @param seed Random seed.
#' @return A `hixscore_train` list: `model` (for [hixscore()]), `E_baseline`, `C_baseline`, `M_model`, `programs`,
#'   `scans` (the M rank scan and the BDM stability table), `intermediate` (reusable BDM results), `status`
#'   (state, time and error of each step) and `training` (sample-level training scores).
#' @export
hixscore_train <- function(epi = NULL, tme = NULL, labels = NULL, steps = "ALL",
                           levels = c("NT", "Inflammatory", "Malignant"),
                           M_model = "default", programs = "default", M_k = "auto", bdm_k = NULL, bdm_drop = integer(),
                           bdm_state_col = "MP_state", bdm_args = list(), bdm_maxK = 10, bdm_reps = 200,
                           M_ranks = 2:10, previous = NULL, on_error = c("continue", "stop"),
                           sample_col = "sample_id", tissue_col = "tissue", series_col = "series",
                           nt_label = "NT", state_col = "cell_state", min_cells_epi = 30, min_cells_tme = 200,
                           add_mp_scores = TRUE, seed = 123, mp_list = NULL, M_W = NULL) {
  .check_seurat()
  on_error <- match.arg(on_error); steps_in <- steps; steps <- .norm_steps(steps); stopped_for_k <- FALSE
  if (!is.null(previous) && !inherits(previous, "hixscore_train")) stop("`previous` must be the result of hixscore_train().", call. = FALSE)
  mp_cols <- NULL
  if (!is.null(mp_list)) { mp_list <- validate_MP_list(mp_list); mp_cols <- names(mp_list) }
  mp_use <- if (is.null(mp_list)) default_MP_list() else mp_list
  if (!is.null(M_W) && !(identical(M_model, "default") || identical(M_model, "derive")))
    stop("Give either `M_W` or an M-layer model object in `M_model`, not both.", call. = FALSE)
  derive_M <- identical(M_model, "derive") || !is.null(M_W); derive_B <- identical(programs, "derive")
  # steps that only make sense when the component is derived from the data
  if (!derive_M) steps <- setdiff(steps, c("M_scan", "M"))
  if (!is.null(M_W)) steps <- setdiff(steps, "M_scan")
  if (!derive_B) steps <- setdiff(steps, c("BDM_scan", "BDM"))
  # "ALL" with `previous` runs only the steps that are not finished yet; named steps are always rerun
  if (any(toupper(as.character(steps_in)) == "ALL") && !is.null(previous$status))
    steps <- setdiff(steps, previous$status$layer[previous$status$status %in% c("done", "reused")])
  # The number of BDM modules cannot be chosen automatically: it is read from the stability table of the scan.
  # Requested explicitly ("BDM"), a missing bdm_k is an error before anything is computed; under "ALL" the run
  # stops after the scan and says what to do next.
  all_requested <- any(toupper(as.character(steps_in)) == "ALL")
  have_scan <- !is.null(previous$scans$bdm)
  if ("BDM" %in% steps && is.null(bdm_k)) {
    if (all_requested) { steps <- setdiff(steps, c("BDM", "C", "fusion")); stopped_for_k <- TRUE }
    else stop("`bdm_k` (number of BDM modules) is required for the \"BDM\" step. ",
              if (have_scan) "Choose it from `previous$scans$bdm$stability`." else
                "Run steps = \"BDM_scan\" first, inspect `$scans$bdm$stability`, then pass `bdm_k` and `previous`.", call. = FALSE)
  }
  if ("fusion" %in% steps) {
    if (is.null(labels)) stop("The \"fusion\" step needs `labels` (data.frame with sample_id and label).", call. = FALSE)
    if (!all(c("sample_id", "label") %in% names(labels))) stop("`labels` needs columns sample_id and label.", call. = FALSE)
    if (!nt_label %in% levels) stop("`levels` must contain the NT label '", nt_label, "'.", call. = FALSE)
  }

  prev <- if (is.null(previous)) list() else unclass(previous)
  R <- list(M_model = prev$M_model, programs = prev$programs, E_baseline = prev$E_baseline, C_baseline = prev$C_baseline,
            thr = prev$thr, model = prev$model, training = prev$training)
  scans <- if (is.null(prev$scans)) list() else prev$scans
  inter <- if (is.null(prev$intermediate)) list() else prev$intermediate
  status <- stats::setNames(lapply(.TRAIN_STEPS, function(s) .status_row(s, "not run")), .TRAIN_STEPS)
  for (s in setdiff(.TRAIN_STEPS, steps)) if (!is.null(prev$status) && prev$status$status[prev$status$layer == s] %in% c("done", "reused"))
    status[[s]] <- .status_row(s, "reused")
  if (is.list(M_model) && !derive_M) R$M_model <- M_model
  if (identical(M_model, "default") && is.null(R$M_model)) R$M_model <- default_M_model()
  if (is.list(programs) && !derive_B) R$programs <- validate_BDM_programs(programs)
  if (identical(programs, "default") && is.null(R$programs)) R$programs <- default_BDM_programs()

  go <- function(step, expr) {
    r <- .run_layer(step, expr, on_error); status[[step]] <<- .status_row(step, r$status, r$seconds, r$error); r }
  need_epi <- function(s) if (is.null(epi)) stop("Step ", s, " needs `epi`.", call. = FALSE)
  need_tme <- function(s) if (is.null(tme)) stop("Step ", s, " needs `tme`.", call. = FALSE)
  md <- if (!is.null(epi)) .get_meta(epi) else NULL
  failed_dep <- function(...) { d <- c(...); d[vapply(d, function(x) status[[x]]$status %in% c("failed", "skipped"), logical(1))] }
  skip <- function(step, why) { status[[step]] <<- .status_row(step, "skipped", NA_real_, why)
                                if (on_error == "stop") stop("Step ", step, " skipped: ", why, call. = FALSE) }

  # ---- M layer ----
  M_composition <- function() { need_tme("M"); build_M_composition(tme, sample_col, state_col, tissue_col = if (!is.null(tme) && tissue_col %in% names(.get_meta(tme))) tissue_col else NULL, min_cells = min_cells_tme) }
  if ("M_scan" %in% steps) { r <- go("M_scan", scan_M_rank(M_composition(), ranks = M_ranks, seed = seed)); if (r$status == "done") scans$M <- r$value }
  if ("M" %in% steps) {
    r <- go("M", {
      P <- M_composition(); k <- if (is.null(M_W)) M_k else NULL
      if (identical(k, "auto")) { k <- attr(if (!is.null(scans$M)) scans$M else scan_M_rank(P, ranks = M_ranks, seed = seed), "optimal_rank")
        if (is.na(k)) stop("Automatic rank selection found no rank (see scan_M_rank()); set `M_k` manually.") }
      if (!is.null(M_W)) build_M_ecotypes(P, W = M_W, seed = seed) else build_M_ecotypes(P, k = k, seed = seed) })
    if (r$status == "done") R$M_model <- r$value
  }

  # ---- BDM gene programs ----
  if ("BDM_scan" %in% steps) {
    r <- go("BDM_scan", {
      need_epi("BDM_scan")
      if (!bdm_state_col %in% names(md)) stop("`epi` lacks the epithelial state column '", bdm_state_col, "'. Annotate with annotate_cells(minor_for = \"all\") and pass bdm_state_col = \"predicted.minor\".")
      common <- list(sample_col = sample_col, state_col = bdm_state_col, tissue_col = tissue_col, series_col = series_col)
      # bdm_args may hold arguments of the three BDM functions and of bdm_pseudobulk(); an unknown name is an error
      bdm_pb <- c("min_cells", "min_logcpm", "exclude_patterns", "subsite_col", "assay")
      bdm_allowed <- list(cg = c("train_cohorts", "max_var", "min_freq", "min_nt_pseudobulk", "workers", bdm_pb),
                          dg = c("train_cohorts", "min_fc", "min_consistency", "min_disease_types", "min_freq", bdm_pb),
                          fd = c("train_cohorts", "max_na", bdm_pb))
      bad <- setdiff(names(bdm_args), unique(unlist(bdm_allowed)))
      if (length(bad)) stop("Unknown `bdm_args`: ", paste(bad, collapse = ", "), ". See ?bdm_constrained_genes, ?bdm_deviation_genes, ?bdm_full_deviation, ?bdm_pseudobulk.")
      pick_args <- function(which) bdm_args[intersect(names(bdm_args), bdm_allowed[[which]])]
      cg <- if (!is.null(inter$constrained)) inter$constrained else do.call(bdm_constrained_genes, c(list(epi, nt_label = nt_label), common, pick_args("cg")))
      dg <- if (!is.null(inter$deviation)) inter$deviation else do.call(bdm_deviation_genes, c(list(epi, cg$genes, nt_label = nt_label), common, pick_args("dg")))
      dev <- if (!is.null(inter$dev)) inter$dev else do.call(bdm_full_deviation, c(list(epi, dg$genes, nt_label = nt_label), common, pick_args("fd")))
      sc <- bdm_scan_modules(dev, maxK = bdm_maxK, reps = bdm_reps, seed = seed)
      list(constrained = cg, deviation = dg, dev = dev, scan = sc) })
    if (r$status == "done") { inter <- utils::modifyList(inter, r$value[c("constrained", "deviation", "dev", "scan")]); scans$bdm <- r$value$scan }
  }
  if ("BDM" %in% steps) {
    if (is.null(inter$dev) || is.null(inter$deviation)) skip("BDM", "run the \"BDM_scan\" step first (no deviation matrix available)")
    else { r <- go("BDM", build_BDM_geneset(inter$dev, inter$deviation$direction, k = bdm_k, drop_modules = bdm_drop, seed = seed, scan = inter$scan))
           if (r$status == "done") R$programs <- r$value }
  }

  # ---- E and C baselines ----
  if ("E" %in% steps) {
    r <- go("E", { need_epi("E"); e <- if (add_mp_scores) compute_MP_scores(epi, mp_list = mp_use) else epi
                   fit_E_baseline(.get_meta(e), sample_col, tissue_col, nt_label, min_cells = min_cells_epi, mp_cols = mp_cols) })
    if (r$status == "done") R$E_baseline <- r$value
  }
  if ("C" %in% steps) {
    if (derive_B && is.null(R$programs)) skip("C", "BDM programs are not available (run the BDM steps or supply `programs`)")
    else {
      r <- go("C", {
        need_epi("C"); n <- table(md[[sample_col]]); e <- epi[, md[[sample_col]] %in% names(n)[n >= min_cells_epi]]
        cb <- fit_C_baseline(e, R$programs, tissue_col, nt_label)
        e <- compute_scHIX(e, baseline = cb, tissue_col = tissue_col, nt_label = nt_label)
        list(baseline = cb, thr = nt_scHIX_thresholds(e, tissue_col = tissue_col, nt_label = nt_label), epi = e) })
      if (r$status == "done") { R$C_baseline <- r$value$baseline; R$thr <- r$value$thr; inter$epi_scHIX <- r$value$epi }
    }
  }

  # ---- fusion model ----
  if ("fusion" %in% steps) {
    dep <- failed_dep("E", "C", if (derive_M) "M")
    miss <- c(E = is.null(R$E_baseline), C = is.null(R$C_baseline) || is.null(R$thr), M = is.null(R$M_model))
    if (length(dep)) skip("fusion", paste0("depends on failed/skipped step(s): ", paste(dep, collapse = ", ")))
    else if (any(miss)) skip("fusion", paste0("missing component(s): ", paste(names(miss)[miss], collapse = ", "), " (run those steps first or pass `previous`)"))
    else {
      r <- go("fusion", {
        need_epi("fusion"); need_tme("fusion")
        n <- table(md[[sample_col]]); keep <- md[[sample_col]] %in% names(n)[n >= min_cells_epi]
        e <- if (!is.null(inter$epi_scHIX)) inter$epi_scHIX else compute_scHIX(epi[, keep], baseline = R$C_baseline, tissue_col = tissue_col, nt_label = nt_label)
        e0 <- if (add_mp_scores) compute_MP_scores(epi[, keep], mp_list = mp_use) else epi[, keep]
        E <- compute_E_seHIX(.get_meta(e0), baseline = R$E_baseline, sample_col = sample_col, tissue_col = tissue_col, min_cells = min_cells_epi)
        Cf <- .sample_C_features(.get_meta(e), "scHIX_cosine", sample_col, R$thr)
        M <- compute_M_CM_scores(tme, R$M_model, sample_col, state_col, min_cells_tme)
        lab <- labels; if (series_col %in% names(md) && !"series" %in% names(lab)) lab$series <- md[[series_col]][match(lab$sample_id, md[[sample_col]])]
        tab <- build_training_table(E, Cf, M, lab)
        message("Training on ", nrow(tab), " samples with E, C and M scores: ", paste(names(table(tab$label)), table(tab$label), collapse = ", "))
        fit_sysHIX_model(tab, levels, R$thr, seed = seed) })
      if (r$status == "done") { R$model <- r$value; R$training <- attr(r$value, "training") }
    }
  }

  if (stopped_for_k) message("BDM_scan finished. Choose the number of modules `bdm_k` from `$scans$bdm$stability`, then rerun with ",
                             "steps = \"ALL\", bdm_k = <k> and previous = <this result>; the finished steps are reused.")
  structure(list(model = R$model, E_baseline = R$E_baseline, C_baseline = R$C_baseline, M_model = R$M_model,
                 programs = R$programs, mp_list = if (!is.null(mp_list)) mp_list else prev$mp_list, thr = R$thr, scans = scans, intermediate = inter,
                 status = do.call(rbind, unname(status[.TRAIN_STEPS])), training = R$training),
            class = "hixscore_train")
}

#' @export
print.hixscore_train <- function(x, ...) {
  cat("<hixscore_train>\n"); print(x$status, row.names = FALSE); invisible(x)
}
