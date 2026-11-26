# Independent oral validation cohort: download the two rds files and check hixscore(mode = "pretrained")

.validation_registry <- list(
  counts = list(file = "validation_counts.rds",
                url = "https://comics.med.sustech.edu.cn/track/scRNA-ref/validation_counts.rds",
                md5 = "3de28ed338875f7b18d94037b7167dab", size = 115671904),
  meta = list(file = "validation_meta.rds",
              url = "https://comics.med.sustech.edu.cn/track/scRNA-ref/validation_meta.rds",
              md5 = "c624a1f191989cb478f3186e2491ad01", size = 747516)
)

# Values written by sysHIX_run_Oral_public_new3.R (baseline_mode = "oral") and reproduced by hixscore(mode = "pretrained");
# NT against the disease groups with more than 3 samples, AUC direction fixed so that a higher stage has a higher score
.validation_expected <- data.frame(
  layer = c("E_seHIX_z", "S_C_z", "S_M_z", "sysHIX"),
  n_samples = c(44L, 44L, 49L, 44L),
  n_NT = c(13L, 13L, 13L, 13L),
  n_disease = c(31L, 31L, 36L, 31L),
  AUC = c(0.803970223325062, 0.88833746898263, 0.942307692307692, 0.935483870967742),
  stringsAsFactors = FALSE
)

.validation_cache_dir <- function() {
  d <- tools::R_user_dir("HIXscore", which = "cache")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

# AUC of class 1 scoring higher, identical to pROC::auc(roc(y, score, levels = c(0, 1), direction = "<"))
.auc_rank <- function(y, score) {
  ok <- !is.na(y) & !is.na(score)
  y <- y[ok]; score <- score[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(score)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# the two paths given as a named list/vector (`counts`, `meta`), as returned by download_validation_data()
.validation_paths <- function(x) {
  if (!is.null(names(x)) && all(c("counts", "meta") %in% names(x))) return(list(counts = x[["counts"]], meta = x[["meta"]]))
  NULL
}

.read_if_path <- function(x) if (is.character(x) && length(x) == 1) readRDS(x) else x

#' Download the validation data
#'
#' Downloads the two files of the independent oral validation cohort into the user cache
#' (`tools::R_user_dir("HIXscore", "cache")`) and verifies their md5 sums:
#'
#' * `validation_counts.rds`: raw counts (`dgCMatrix`, genes x epithelial cells) of the epithelial cells;
#' * `validation_meta.rds`: one row per cell with the anonymous sample id (`S01`...), the tissue
#'   (`NT`, `OSF`, `OLP`, `CHC`, `OLK`, `Tumor`), the major cell type and the cell state.
#'
#' Cells are numbered consecutively; no barcode, cohort or patient information is kept. Cached files with the correct
#' md5 sum are reused. The md5 sum detects a damaged or changed download; it does not establish who produced the file.
#'
#' @param dest Directory for the files (default: the HIXscore cache).
#' @param force Download again even if verified copies are cached.
#' @param timeout Download time limit in seconds.
#' @return Named character vector with the paths `counts` and `meta`.
#' @export
download_validation_data <- function(dest = NULL, force = FALSE, timeout = 3600) {
  if (is.null(dest)) dest <- .validation_cache_dir()
  if (!dir.exists(dest)) dir.create(dest, recursive = TRUE)
  old <- options(timeout = max(timeout, getOption("timeout"))); on.exit(options(old))
  out <- character(0)
  for (nm in names(.validation_registry)) {
    r <- .validation_registry[[nm]]; f <- file.path(dest, r$file)
    if (file.exists(f) && !force && unname(tools::md5sum(f)) == r$md5) { out[[nm]] <- f; next }
    tmp <- paste0(f, ".part"); on.exit(unlink(tmp), add = TRUE)
    utils::download.file(r$url, tmp, mode = "wb", quiet = FALSE)
    if (unname(tools::md5sum(tmp)) != r$md5)
      stop("md5 mismatch for the downloaded file '", r$file, "'; file removed. The file on the server differs from the ",
           "one the expected scores belong to.", call. = FALSE)
    file.rename(tmp, f)
    out[[nm]] <- f
  }
  out
}

#' Validate a HIXscore installation on the independent oral cohort
#'
#' Runs [hixscore()] with `mode = "pretrained"` on the validation data ([download_validation_data()]) and compares the
#' AUCs of the four layer scores with the values of the reference pipeline
#' (`sysHIX_run_Oral_public_new3.R`, oral baseline mode). The AUC is NT against the disease groups with more than 3
#' samples (`OLP`, `OSF`, `Tumor`), with a higher score for the higher stage.
#'
#' A correct installation reproduces the reference AUCs (E 0.804, C 0.888, M 0.942 and sysHIX 0.935 with 13 NT samples)
#' and the same number of NT and disease samples per layer.
#'
#' @param counts Path to `validation_counts.rds`, the matrix itself, or a list/vector of the two paths returned by
#'   [download_validation_data()]. `NULL` downloads the data.
#' @param meta Path to `validation_meta.rds`, the data frame itself, or the second path of the validation data.
#' @param tolerance Accepted absolute difference between the observed and the reference AUC.
#' @param verbose Print the comparison table.
#' @param ... Passed to [hixscore()].
#' @return A data.frame with one row per layer: `layer`, the expected and observed sample numbers (`n_samples`,
#'   `n_NT`, `n_disease`), `AUC_expected`, `AUC` and `pass`; attribute `"all_pass"` is `TRUE` when every layer matches.
#' @export
validate_hixscore <- function(counts = NULL, meta = NULL, tolerance = 0.005, verbose = TRUE, ...) {
  .check_seurat()
  if (is.null(counts) && is.null(meta)) {
    p <- download_validation_data(); counts <- p[["counts"]]; meta <- p[["meta"]]
  } else if (is.null(meta)) {
    p <- .validation_paths(counts); if (is.null(p)) stop("Pass both files: `counts` and `meta`.", call. = FALSE)
    counts <- p$counts; meta <- p$meta
  }
  counts <- .read_if_path(counts); meta <- .read_if_path(meta)
  if (!inherits(counts, "Matrix")) stop("`counts` must be a dgCMatrix-shaped count matrix.", call. = FALSE)
  need <- c("sample_id", "tissue", "major", "cell_state")
  if (!is.data.frame(meta) || !all(need %in% names(meta)))
    stop("`meta` must be a data.frame with the columns ", paste(need, collapse = ", "), ".", call. = FALSE)
  cells <- rownames(meta)
  if (!identical(colnames(counts), cells[meta$major == "Epithelial"]))
    stop("The columns of `counts` must be the epithelial cells of `meta`, in the same order (cell numbers).", call. = FALSE)

  epi <- Seurat::CreateSeuratObject(counts = counts, project = "HIXscore validation", min.cells = 0, min.features = 0)
  md <- meta[match(colnames(epi), cells), c("sample_id", "tissue")]
  epi$sample_id <- as.character(md$sample_id); epi$tissue <- as.character(md$tissue)
  epi <- Seurat::NormalizeData(epi, verbose = FALSE)
  tme <- data.frame(sample_id = as.character(meta$sample_id[meta$major != "Epithelial"]),
                    cell_state = as.character(meta$cell_state[meta$major != "Epithelial"]), stringsAsFactors = FALSE)
  if (verbose) message("Running hixscore(mode = 'pretrained') on ", ncol(epi), " epithelial and ", nrow(tme),
                       " immune/stromal cells...")
  res <- suppressWarnings(hixscore(epi, tme, mode = "pretrained", sample_col = "sample_id", tissue_col = "tissue", ...))
  s <- res$samples

  out <- do.call(rbind, lapply(.validation_expected$layer, function(v) {
    d <- data.frame(tissue = as.character(s$tissue), score = s[[v]]); d <- d[!is.na(d$score), ]
    n <- table(d$tissue); keep <- names(n)[n > 3]; d <- d[d$tissue %in% keep, ]
    y <- as.integer(d$tissue != "NT")
    data.frame(layer = v, n_samples = nrow(d), n_NT = sum(y == 0), n_disease = sum(y == 1),
               AUC = .auc_rank(y, d$score), stringsAsFactors = FALSE)
  }))
  out$AUC_expected <- .validation_expected$AUC[match(out$layer, .validation_expected$layer)]
  k <- match(out$layer, .validation_expected$layer)
  out$pass <- abs(out$AUC - out$AUC_expected) <= tolerance & out$n_NT == .validation_expected$n_NT[k] &
    out$n_disease == .validation_expected$n_disease[k]
  out <- out[, c("layer", "n_samples", "n_NT", "n_disease", "AUC_expected", "AUC", "pass")]
  if (verbose) {
    print(format(out, digits = 4), row.names = FALSE)
    if (all(out$pass)) message("OK: the layer scores reproduce the reference AUCs (tolerance ", tolerance, ").")
  }
  if (!all(out$pass))
    warning("Layer(s) outside the expected values: ", paste(out$layer[!out$pass], collapse = ", "),
            ". Check the Seurat version (4.3-4.4), the HIXscore version and the reproducibility of the environment.",
            call. = FALSE)
  attr(out, "all_pass") <- all(out$pass)
  attr(out, "samples") <- s
  out
}
