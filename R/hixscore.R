# hixscore(): run the E, C, M layers and the sysHIX fusion separately or together.
# Every layer is a self-contained step that uses only its own inputs and filters; results are kept in the
# returned object and can be passed back as `previous` so that only the missing or failed layers are rerun.
# Samples are intersected only in the sysHIX step.

.LAYERS <- c("E", "C", "M", "sysHIX")

.norm_layers <- function(layers) {
  layers <- unique(as.character(layers))
  if ("ALL" %in% toupper(layers)) return(.LAYERS)
  layers <- ifelse(tolower(layers) == "syshix", "sysHIX", toupper(layers))
  bad <- setdiff(layers, .LAYERS)
  if (length(bad)) stop("Unknown layer(s): ", paste(bad, collapse = ", "), ". Use any of E, C, M, sysHIX, or \"ALL\".", call. = FALSE)
  .LAYERS[.LAYERS %in% layers]
}

# run one layer, never letting an error escape unless on_error == "stop"
.run_layer <- function(name, expr, on_error) {
  t0 <- Sys.time()
  r <- tryCatch(list(value = expr, error = NA_character_), error = function(e) list(value = NULL, error = conditionMessage(e)))
  if (!is.na(r$error) && on_error == "stop") stop("Layer ", name, " failed: ", r$error, call. = FALSE)
  if (!is.na(r$error)) warning("Layer ", name, " failed: ", r$error, call. = FALSE)
  list(value = r$value, status = if (is.na(r$error)) "done" else "failed", error = r$error,
       seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
}

.status_row <- function(layer, status, seconds = NA_real_, note = NA_character_)
  data.frame(layer = layer, status = status, seconds = seconds, note = note, stringsAsFactors = FALSE)

#' Compute the HIXscore layers
#'
#' Computes the E layer (seHIX), the C layer (scHIX -> S_C), the M layer (CM -> S_M) and the sysHIX
#' fusion, either all at once or one layer at a time. Each layer is scored on the samples that pass its
#' own filters and needs only its own input: E and C need `epi`, M needs `tme` and the number of epithelial
#' cells of each sample (taken from `epi` or given as `n_epi`), and `sysHIX` fuses the three layer results
#' (computed now or taken from `previous`). Samples are intersected only in the `sysHIX` step.
#'
#' A layer that fails does not stop the others (`on_error = "continue"`): the error is recorded in the
#' `status` table and the layers that depend on it are skipped. Pass the returned object as `previous` to
#' rerun only selected layers; layers not requested are taken from `previous`.
#'
#' @param epi Epithelial Seurat object (log-normalised) with the sample column; `tissue_col` is needed
#'   when `mode = "custom_baseline"`.
#' @param tme Non-epithelial cells (Seurat object or data.frame) with the sample column and `state_col`.
#' @param layers Layers to run: any of `"E"`, `"C"`, `"M"`, `"sysHIX"`, or `"ALL"`.
#' @param mode `"pretrained"`: shipped oral baselines and model. `"custom_baseline"`: E and C baselines and the
#'   scHIX thresholds are refitted on the NT samples of `epi`.
#' @param previous Result of an earlier `hixscore()` call whose layer results and components are reused.
#' @param n_epi Optional named numeric vector (names = sample ids) of epithelial cells per sample, used by the
#'   M layer when `epi` is not supplied.
#' @param on_error `"continue"` (default) records a failed layer and goes on; `"stop"` stops at the first error.
#' @param sample_col,tissue_col,state_col Metadata columns for sample id, tissue label and immune/stromal cell state.
#' @param nt_label Label of normal tissue in `tissue_col`.
#' @param min_cells_epi,min_cells_tme Minimum epithelial / immune+stromal cells per sample.
#' @param min_epi_M Minimum epithelial cells of a sample for its M-layer score; the default 0 keeps all samples
#'   (those at or below the training range are flagged in `M_epi_below_training_min`).
#' @param add_mp_scores Compute MP scores on `epi` (set `FALSE` if the MP columns already exist).
#' @param model sysHIX model (default shipped oral model, or the `model` from [hixscore_train()]).
#' @param E_baseline,C_baseline,M_model Components from [hixscore_train()] to use instead of the shipped ones.
#' @param programs Optional user BDM gene sets (named list); requires a model from [hixscore_train()].
#' @param mp_list Optional user E-layer gene sets (named list, see [validate_MP_list()]) used instead of
#'   [default_MP_list()] for the MP scores. In `mode = "custom_baseline"` the E baseline is refitted on these sets. In
#'   `mode = "pretrained"` the shipped baseline and model are bound to the shipped MPs, so the names of `mp_list` must
#'   equal those of [default_MP_list()] (an error otherwise); with an `E_baseline` from [hixscore_train()] the names
#'   must include the MPs of that baseline.
#' @return A `hixscore` list: `samples` (one row per sample seen by any layer with the independent E, S_C and S_M
#'   scores, z-scores, cell counts, `has_*` flags and sysHIX for samples with all three layers), `E`, `C`, `M`
#'   (per-layer tables), `sysHIX` (the intersection with class probabilities), `status` (state, time and error of
#'   each layer), `components` (baselines, thresholds and model used), `epi` (epithelial object with scHIX),
#'   `n_epi`, `info`, `mode`.
#' @export
hixscore <- function(epi = NULL, tme = NULL, layers = "ALL", mode = c("pretrained", "custom_baseline"),
                     previous = NULL, n_epi = NULL, on_error = c("continue", "stop"),
                     sample_col = "sample_id", tissue_col = "tissue", nt_label = "NT", state_col = "cell_state",
                     min_cells_epi = 30, min_cells_tme = 200, min_epi_M = 0, add_mp_scores = TRUE,
                     model = default_sysHIX_model(), E_baseline = NULL, C_baseline = NULL,
                     M_model = default_M_model(), programs = NULL, mp_list = NULL) {
  .check_seurat()
  mode <- match.arg(mode); on_error <- match.arg(on_error); layers <- .norm_layers(layers)
  if (!is.null(previous) && !inherits(previous, "hixscore")) stop("`previous` must be the result of hixscore().", call. = FALSE)
  if (!is.null(epi)) check_input(epi, NULL, sample_col, tissue_col, state_col, mode, nt_label)
  if (!is.null(tme)) {
    tm <- .get_meta(tme)
    if (!all(c(sample_col, state_col) %in% names(tm)))
      stop("Non-epithelial metadata needs columns '", sample_col, "' and '", state_col, "'.", call. = FALSE)
  }
  if (!is.null(programs)) {
    programs <- validate_BDM_programs(programs)
    if (!inherits(model, "sysHIX_model"))
      stop("Custom `programs` need a model trained on them: run hixscore_train() and pass its `model`, `E_baseline`, ",
           "`C_baseline` and `M_model` to hixscore().", call. = FALSE)
  }
  if (mode == "custom_baseline")
    message("custom_baseline: E and C baselines (and scHIX tertile thresholds) are refitted on your NT samples",
            if (inherits(model, "sysHIX_model")) "; the user-trained model is used." else
              "; layer-to-score mappings and the fusion model stay at the oral training values.")

  mp_cols <- NULL
  if (!is.null(mp_list)) {
    mp_list <- validate_MP_list(mp_list); mp_cols <- names(mp_list)
    if (!is.null(E_baseline)) {
      bm <- .as_E_baseline(E_baseline)$mp_cols
      if (!all(bm %in% mp_cols)) stop("`mp_list` lacks MP(s) of the supplied E_baseline: ", paste(setdiff(bm, mp_cols), collapse = ", "), call. = FALSE)
    } else if (mode == "pretrained" && !setequal(mp_cols, names(default_MP_list())))
      stop("In mode 'pretrained' the shipped E baseline and model are bound to the shipped MPs, so `mp_list` must have the names of default_MP_list(). ",
           "Use mode = 'custom_baseline' (E baseline refitted on your NT samples) or a model from hixscore_train().", call. = FALSE)
  }

  prev <- if (is.null(previous)) list() else unclass(previous)
  res <- list(E = prev$E, C = prev$C, M = prev$M, sysHIX = prev$sysHIX, epi = prev$epi)
  comp <- if (is.null(prev$components)) list() else prev$components
  status <- list()
  mark <- function(layer, r) { status[[layer]] <<- .status_row(layer, r$status, r$seconds, r$error) }
  for (l in setdiff(.LAYERS, layers))
    status[[l]] <- if (!is.null(res[[l]])) .status_row(l, "reused") else .status_row(l, "not run")

  md <- if (!is.null(epi)) .get_meta(epi) else NULL
  need_epi <- function(l) if (is.null(epi)) stop("Layer ", l, " needs `epi`.", call. = FALSE)

  # shared setup: the model used by every layer; MP scores are added when the E layer is run
  comp$model <- model
  if ("E" %in% layers && add_mp_scores && !is.null(epi)) { epi <- compute_MP_scores(epi, mp_list = if (is.null(mp_list)) default_MP_list() else mp_list); md <- .get_meta(epi) }

  # ---- E layer ----
  if ("E" %in% layers) {
    r <- .run_layer("E", {
      need_epi("E")
      E_base <- if (!is.null(E_baseline)) E_baseline else if (mode == "pretrained") system.file("extdata", "E_seHIX_baseline_oral.rds", package = "HIXscore") else NULL
      E <- compute_E_seHIX(md, baseline = E_base, sample_col = sample_col, tissue_col = tissue_col, nt_label = nt_label, min_cells = min_cells_epi, mp_cols = mp_cols)
      compute_E_z(E, model)
    }, on_error)
    mark("E", r); if (r$status == "done") res$E <- r$value
  }

  # ---- C layer ----
  if ("C" %in% layers) {
    r <- .run_layer("C", {
      need_epi("C")
      C_base <- if (!is.null(C_baseline)) C_baseline else if (mode == "pretrained") system.file("extdata", "scHIX_baseline_oral.rds", package = "HIXscore") else NULL
      if (!is.null(programs) && is.null(C_baseline)) C_base <- fit_C_baseline(epi, programs, tissue_col, nt_label)
      # samples below the cell threshold are dropped BEFORE module scoring (AddModuleScore picks control genes from the cells it gets)
      n_per <- table(md[[sample_col]]); keep <- md[[sample_col]] %in% names(n_per)[n_per >= min_cells_epi]
      if (!any(keep)) stop("No sample has at least ", min_cells_epi, " epithelial cells.")
      e2 <- compute_scHIX(epi[, keep], baseline = C_base, tissue_col = tissue_col, nt_label = nt_label)
      mdl <- model
      if (mode == "custom_baseline" && is.null(C_baseline)) {
        thr <- nt_scHIX_thresholds(e2, tissue_col = tissue_col, nt_label = nt_label)
        mdl$scHIX_thresholds$thr_low <- thr$thr_low; mdl$scHIX_thresholds$thr_high <- thr$thr_high
      }
      list(C = compute_S_C(e2, model = mdl, sample_col = sample_col, min_cells = min_cells_epi), epi = e2, model = mdl)
    }, on_error)
    mark("C", r)
    if (r$status == "done") { res$C <- r$value$C; res$epi <- r$value$epi; comp$model <- r$value$model }
  }

  # ---- epithelial cells per sample (M layer covariate) ----
  if (!is.null(epi)) n_epi <- table(as.character(md[[sample_col]]))
  if (is.null(n_epi)) n_epi <- prev$n_epi

  # ---- M layer ----
  if ("M" %in% layers) {
    r <- .run_layer("M", {
      if (is.null(tme)) stop("Layer M needs `tme`.", call. = FALSE)
      if (is.null(n_epi)) stop("Layer M needs the epithelial cell number of each sample: supply `epi`, `n_epi` or `previous`.", call. = FALSE)
      M <- compute_M_CM_scores(tme, model = M_model, sample_col = sample_col, state_col = state_col, min_cells = min_cells_tme)
      compute_S_M(M, n_epi, model, min_epi = min_epi_M)
    }, on_error)
    mark("M", r); if (r$status == "done") res$M <- r$value
    if (r$status == "done" && any(r$value$M_epi_below_training_min))
      message(sum(r$value$M_epi_below_training_min), " sample(s) have <= 10 epithelial cells (", sum(r$value$n_epi_cells == 0),
              " without any): their S_M is kept; see `M_epi_below_training_min`.")
  }

  # ---- sysHIX fusion: the only step that intersects samples ----
  if ("sysHIX" %in% layers) {
    miss <- c(E = is.null(res$E), C = is.null(res$C), M = is.null(res$M)); miss <- names(miss)[miss]
    failed <- vapply(status[miss], function(s) s$status == "failed", logical(1))
    if (length(miss)) {
      status$sysHIX <- .status_row("sysHIX", "skipped", NA_real_, paste0("missing layer(s): ", paste(miss, collapse = ", "),
                                    if (any(failed)) paste0(" (failed: ", paste(miss[failed], collapse = ", "), ")") else ""))
      if (on_error == "stop") stop("sysHIX needs the layers ", paste(miss, collapse = ", "), ".", call. = FALSE)
    } else {
      r <- .run_layer("sysHIX", predict_sysHIX(res$E, res$C, res$M, comp$model), on_error)
      mark("sysHIX", r); if (r$status == "done") res$sysHIX <- r$value
    }
  }

  tis_now <- c(if (!is.null(md) && tissue_col %in% names(md)) .first_per_sample(md[[sample_col]], md[[tissue_col]]),
               if (!is.null(tme) && tissue_col %in% names(.get_meta(tme))) .first_per_sample(.get_meta(tme)[[sample_col]], .get_meta(tme)[[tissue_col]]))
  tissue_map <- c(tis_now, prev$tissue_map); tissue_map <- tissue_map[!duplicated(names(tissue_map))]
  samples <- .assemble_samples(res, tissue_map, n_epi)
  info <- list(mode = mode, min_cells_epi = min_cells_epi, min_cells_tme = min_cells_tme,
               n_nt_E_baseline = { n <- if (!is.null(res$E)) attr(res$E, "n_nt_baseline") else NULL; if (is.null(n)) NA_integer_ else n },
               n_M_epi_below_training_min = if (!is.null(res$M)) sum(res$M$M_epi_below_training_min) else NA_integer_)
  structure(list(samples = samples, E = res$E, C = res$C, M = res$M, sysHIX = res$sysHIX,
                 status = do.call(rbind, unname(status[.LAYERS])), components = comp, epi = res$epi,
                 n_epi = n_epi, tissue_map = tissue_map, info = info, mode = mode), class = "hixscore")
}

# one row per sample seen by any layer
.first_per_sample <- function(sample, value) {
  v <- as.character(value); names(v) <- as.character(sample); v[!duplicated(names(v))]
}

.assemble_samples <- function(res, tissue_map, n_epi) {
  E <- res$E; C <- res$C; M <- res$M; S <- res$sysHIX
  ids <- Reduce(union, list(E$sample_id, C$sample_id, M$sample_id))
  if (is.null(ids)) ids <- character()
  out <- data.frame(sample_id = ids, stringsAsFactors = FALSE)
  if (length(tissue_map)) out$tissue <- unname(tissue_map[ids])
  na_col <- function(type) rep(type, length(ids))
  mt <- function(tab, col, type) if (is.null(tab)) na_col(type) else tab[[col]][match(ids, tab$sample_id)]
  ne <- if (is.null(n_epi)) na_col(NA_integer_) else as.integer(n_epi[ids])
  if (!is.null(n_epi)) ne[is.na(ne)] <- 0L
  out$n_epi_cells <- ne
  out$n_tme_cells <- as.integer(mt(M, "n_cells", NA_integer_))
  out$has_E <- ids %in% E$sample_id; out$has_C <- ids %in% C$sample_id; out$has_M <- ids %in% M$sample_id
  out$has_sysHIX <- ids %in% S$sample_id
  out$E_seHIX <- as.numeric(mt(E, "E_seHIX", NA_real_)); out$E_seHIX_z <- as.numeric(mt(E, "E_seHIX_z", NA_real_))
  out$S_C <- as.numeric(mt(C, "S_C", NA_real_));         out$S_C_z <- as.numeric(mt(C, "S_C_z", NA_real_))
  out$S_M <- as.numeric(mt(M, "S_M", NA_real_));         out$S_M_z <- as.numeric(mt(M, "S_M_z", NA_real_))
  out$M_epi_below_training_min <- as.logical(mt(M, "M_epi_below_training_min", NA))
  out$sysHIX <- as.numeric(mt(S, "sysHIX", NA_real_))
  out
}

#' @export
print.hixscore <- function(x, ...) {
  cat("<hixscore> mode:", x$mode, "| samples:", nrow(x$samples), "\n"); print(x$status, row.names = FALSE)
  invisible(x)
}
