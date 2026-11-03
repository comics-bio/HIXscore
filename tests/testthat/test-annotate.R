# Flow and failure paths of annotate_cells() on a tiny synthetic reference (accuracy is NOT tested here)
synth_reference <- function(seed = 21) {
  skip_if_not_installed("Seurat"); set.seed(seed)
  genes <- paste0("G", 1:1500)
  mk <- function(maj, minors, n_per, shift) {
    do.call(cbind, lapply(seq_along(minors), function(j) {
      mu <- rexp(length(genes), 1); mu[shift + (j - 1) * 20 + 1:20] <- mu[shift + (j - 1) * 20 + 1:20] + 6
      matrix(rpois(length(genes) * n_per, rep(mu, n_per)), length(genes), n_per,
             dimnames = list(genes, paste0(maj, "_", minors[j], "_", seq_len(n_per))))
    }))
  }
  A <- mk("A", c("a1", "a2"), 80, 0); B <- mk("B", c("b1", "b2", "b3"), 70, 200)
  m <- cbind(A, B)
  o <- Seurat::CreateSeuratObject(m)
  o$celltype_major <- sub("_.*", "", colnames(m)); o$cell_type2 <- sub("^[AB]_([a-z0-9]+)_.*", "\\1", colnames(m))
  o <- Seurat::NormalizeData(o, verbose = FALSE); o <- Seurat::FindVariableFeatures(o, nfeatures = 600, verbose = FALSE)
  o <- Seurat::ScaleData(o, verbose = FALSE); o <- Seurat::RunPCA(o, npcs = 20, verbose = FALSE)
  o[["harmony"]] <- Seurat::CreateDimReducObject(embeddings = Seurat::Embeddings(o, "pca"), loadings = Seurat::Loadings(o, "pca"),
                                                 key = "harmony_", assay = "RNA")
  list(ref = o, genes = genes, mk = mk)
}

test_that("annotate_cells() returns labels, handles tiny samples, and refuses Seurat >= 5-style misuse", {
  skip_on_cran(); s <- synth_reference()
  ref <- suppressWarnings(prepare_reference(s$ref, major_col = "celltype_major", minor_col = "cell_type2", min_cells_per_type = 20))
  expect_s3_class(ref, "HIXscore_reference")
  expect_true(all(c("A", "B") %in% names(ref$minor)))
  que <- cbind(s$mk("A", c("a1", "a2"), 60, 0), s$mk("B", c("b1", "b2", "b3"), 50, 200)); colnames(que) <- paste0("q", seq_len(ncol(que)))
  q <- Seurat::CreateSeuratObject(que); q$sample_id <- rep(c("big", "tiny"), c(ncol(que) - 12, 12))
  out <- suppressWarnings(suppressMessages(annotate_cells(q, ref, sample_col = "sample_id", min_cells_per_sample = 50,
                                                          minor_for = "all", dims = 1:10, verbose = FALSE)))
  md <- out@meta.data
  expect_true(all(c("predicted.major", "prediction.score.major", "predicted.minor", "annotation.confidence") %in% names(md)))
  expect_true(all(is.na(md$predicted.major[md$sample_id == "tiny"])))               # below min_cells_per_sample: not annotated
  big <- md[md$sample_id == "big", ]; expect_true(mean(!is.na(big$predicted.major)) > 0.99)
  expect_true(all(big$predicted.major %in% c("A", "B")))
  s2 <- annotation_summary(out); expect_true(s2$annotated_fraction < 1 && s2$annotated_fraction > 0.8)
})

test_that("annotate_cells() input checks", {
  skip_on_cran(); s <- synth_reference()
  ref <- suppressWarnings(prepare_reference(s$ref, major_col = "celltype_major", minor_col = "cell_type2", min_cells_per_type = 20))
  q <- Seurat::CreateSeuratObject(matrix(rpois(1500 * 60, 2), 1500, 60, dimnames = list(s$genes, paste0("q", 1:60))))
  expect_error(annotate_cells(q, ref, sample_col = "nope"))
  expect_error(annotate_cells(q, list()), "HIXscore_reference|inherits")
})

test_that("reference helpers", {
  expect_error(load_reference(path = tempfile()), "not found")
  f <- tempfile(fileext = ".rds"); saveRDS(list(a = 1), f); expect_error(load_reference(path = f), "must be a Seurat object")
  expect_true(is.list(default_ref_features()) && "Epithelial" %in% names(default_ref_features()))
})

test_that("load_reference() validates the reference object and loads a prepared one unchanged", {
  skip_on_cran(); s <- synth_reference(); f <- tempfile(fileext = ".rds")
  base <- s$ref; base$celltype_major <- ifelse(base$celltype_major == "A", "Epithelial", base$celltype_major)
  r1 <- base; r1$celltype_major <- NULL; saveRDS(r1, f); expect_error(load_reference(path = f), "lacks metadata column")
  r2 <- base; r2$celltype_major[1] <- NA; saveRDS(r2, f); expect_error(load_reference(path = f), "missing values")
  r3 <- base; r3$celltype_major <- "A"; saveRDS(r3, f); expect_error(load_reference(path = f), "Epithelial")
  r4 <- base; r4@reductions$pca <- NULL; saveRDS(r4, f); expect_error(load_reference(path = f), "no 'pca' reduction")
  prep <- suppressWarnings(prepare_reference(base, major_col = "celltype_major", minor_col = "cell_type2", min_cells_per_type = 20))
  sv <- save_reference(prep, f); expect_true(nzchar(attr(sv, "md5")))
  back <- suppressWarnings(load_reference(path = f)); expect_s3_class(back, "HIXscore_reference"); expect_equal(back$reduction, "pca")
})

test_that("download_reference() rejects a file with the wrong md5 and reuses a verified cache", {
  skip_on_cran()
  src <- tempfile(fileext = ".rds"); saveRDS(1:10, src)
  withr_dir <- tempfile(); dir.create(withr_dir)
  testthat::local_mocked_bindings(.reference_cache_dir = function() withr_dir, .package = "HIXscore")
  reg <- HIXscore:::.reference_registry
  testthat::local_mocked_bindings(.reference_registry = list(HIXscore = list(file = "r.rds", url = paste0("file://", src), md5 = "0000", size = 0)), .package = "HIXscore")
  expect_error(suppressWarnings(download_reference()), "md5 mismatch")
  expect_false(file.exists(file.path(withr_dir, "r.rds")))                      # nothing left behind
  good <- unname(tools::md5sum(src))
  testthat::local_mocked_bindings(.reference_registry = list(HIXscore = list(file = "r.rds", url = paste0("file://", src), md5 = good, size = 0)), .package = "HIXscore")
  p <- suppressMessages(download_reference()); expect_true(file.exists(p)); expect_equal(unname(tools::md5sum(p)), good)
  unlink(src); expect_equal(download_reference(), p)                            # cached copy is reused, the source is not needed
})

test_that("reduction argument is validated", {
  skip_on_cran(); s <- synth_reference()
  ref <- suppressWarnings(prepare_reference(s$ref, major_col = "celltype_major", minor_col = "cell_type2", min_cells_per_type = 20))
  q <- Seurat::CreateSeuratObject(matrix(rpois(1500 * 60, 2), 1500, 60, dimnames = list(s$genes, paste0("q", 1:60))))
  q$sample_id <- "s"
  expect_error(annotate_cells(q, ref, reduction = "umap"), "no reduction 'umap'|Reference has no reduction")
})

test_that("prepare_input() splits epithelial and non-epithelial cells and checks the vocabulary", {
  skip_on_cran(); set.seed(31)
  st <- rownames(default_M_model()$W); genes <- paste0("G", 1:1500); n <- 400
  o <- Seurat::CreateSeuratObject(matrix(rpois(1500 * n, 3), 1500, n, dimnames = list(genes, paste0("c", 1:n))))
  o$sample_id <- rep(c("a", "b"), each = n / 2)
  o$predicted.major <- rep(c("Epithelial", "CD4_T"), length.out = n); o$predicted.minor <- ifelse(o$predicted.major == "CD4_T", st[1], "MP_1_Keratinization")
  p <- suppressMessages(prepare_input(o))
  expect_equal(ncol(p$epi), n / 2); expect_equal(nrow(p$tme), n / 2); expect_true(all(p$tme$cell_state == st[1]))
  expect_equal(p$summary$in_model_state_fraction, 1)
  o$predicted.minor[o$predicted.major == "CD4_T"] <- "not_a_state"
  expect_warning(suppressMessages(prepare_input(o)), "cell state used by the M-layer model")
  o$predicted.major <- NULL
  expect_error(prepare_input(o), "Missing metadata")
  expect_error(prepare_input(o, annotate = TRUE), "needs a `reference`")
})


test_that("reduction = 'pca' prepares a within-type PCA and cannot be requested on a reference prepared with another reduction", {
  skip_on_cran(); s <- synth_reference()
  rp <- suppressWarnings(prepare_reference(s$ref, major_col = "celltype_major", minor_col = "cell_type2",
                                           min_cells_per_type = 20, reduction = "pca"))
  rh <- suppressWarnings(prepare_reference(s$ref, major_col = "celltype_major", minor_col = "cell_type2",
                                           min_cells_per_type = 20, reduction = "harmony"))
  sub_p <- rp$minor[["A"]]
  expect_true("pca" %in% names(sub_p@reductions)); expect_equal(ncol(Seurat::Embeddings(sub_p, "pca")), max(min(50, floor(ncol(sub_p) * 0.5), nrow(sub_p) - 1), 10))
  expect_equal(nrow(Seurat::Embeddings(sub_p, "pca")), ncol(sub_p))                                          # recomputed on the subset
  q <- Seurat::CreateSeuratObject(matrix(rpois(1500 * 60, 2), 1500, 60, dimnames = list(s$genes, paste0("q", 1:60)))); q$sample_id <- "s"
  expect_error(annotate_cells(q, rh, reduction = "pca"), "Re-run prepare_reference")
})

test_that("a failing sample does not stop annotate_cells(), and checkpoints are reused", {
  skip_on_cran(); s <- synth_reference()
  ref <- suppressWarnings(prepare_reference(s$ref, major_col = "celltype_major", minor_col = "cell_type2", min_cells_per_type = 20))
  que <- cbind(s$mk("A", c("a1", "a2"), 40, 0), s$mk("B", c("b1", "b2", "b3"), 30, 200)); colnames(que) <- paste0("q", seq_len(ncol(que)))
  q <- Seurat::CreateSeuratObject(que); q$sample_id <- rep(c("ok1", "bad", "ok2"), length.out = ncol(que))
  q$sample_id[seq_len(60)] <- "ok1"; q$sample_id[61:120] <- "bad"; q$sample_id[121:ncol(q)] <- "ok2"
  ck <- tempfile(); calls <- 0
  real <- HIXscore:::.transfer_major
  testthat::local_mocked_bindings(.transfer_major = function(ref, que, ...) { calls <<- calls + 1
    if (grepl("^q61$", colnames(que)[1])) stop("simulated failure"); real(ref, que, ...) }, .package = "HIXscore")
  out <- suppressWarnings(annotate_cells(q, ref, sample_col = "sample_id", dims = 1:10, checkpoint_dir = ck, verbose = FALSE))
  md <- out@meta.data
  expect_true(all(is.na(md$predicted.major[md$sample_id == "bad"])))             # failed sample keeps NA
  expect_true(mean(!is.na(md$predicted.major[md$sample_id != "bad"])) > 0.99)    # the others are annotated
  f <- attr(out, "annotation_failed"); expect_equal(f$sample, "bad"); expect_match(f$reason, "simulated failure")
  expect_equal(length(list.files(ck)), 2)                                        # only the finished samples are checkpointed
  calls <- 0
  out2 <- suppressWarnings(annotate_cells(q, ref, sample_col = "sample_id", dims = 1:10, checkpoint_dir = ck, verbose = FALSE))
  expect_equal(calls, 1)                                                         # only the failed sample is attempted again
  expect_equal(out2$predicted.major[q$sample_id != "bad"], out$predicted.major[q$sample_id != "bad"])
})
