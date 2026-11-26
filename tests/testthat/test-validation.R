test_that(".auc_rank equals the pROC AUC with direction '<'", {
  skip_if_not_installed("pROC")
  set.seed(1); y <- rep(0:1, c(15, 25)); s <- c(rnorm(15), rnorm(25, 1)); s[3] <- s[20]            # include a tie
  expect_equal(HIXscore:::.auc_rank(y, s), as.numeric(pROC::auc(pROC::roc(y, s, levels = c(0, 1), direction = "<", quiet = TRUE))))
  expect_true(is.na(HIXscore:::.auc_rank(rep(1, 5), 1:5)))
})

test_that("the expected values and the registry are consistent", {
  e <- HIXscore:::.validation_expected
  expect_equal(e$layer, c("E_seHIX_z", "S_C_z", "S_M_z", "sysHIX"))
  expect_true(all(e$AUC > 0.5 & e$AUC < 1))
  r <- HIXscore:::.validation_registry
  expect_setequal(names(r), c("counts", "meta"))
  expect_true(all(vapply(r, function(x) grepl("^[0-9a-f]{32}$", x$md5) && grepl("^https://", x$url), logical(1))))
})

test_that("validate_hixscore() checks its input", {
  skip_if_not_installed("Seurat")
  m <- Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(3, 2), dimnames = list(letters[1:3], c("1", "2")))
  meta <- data.frame(sample_id = c("S01", "S01"), tissue = "NT", major = c("Epithelial", "T"), cell_state = c("a", "b"),
                     row.names = c("1", "2"))
  expect_error(validate_hixscore(m, meta[, c("sample_id", "tissue")]), "columns")
  expect_error(validate_hixscore(m, meta), "epithelial cells")                                      # columns of counts != epithelial cells
  expect_error(validate_hixscore(as.matrix(m), meta), "count matrix")
  expect_error(validate_hixscore(list(a = 1)), "both files")
})

test_that("validate_hixscore() compares the AUCs and warns when they differ", {
  skip_on_cran(); d <- synth_data()
  epi <- d$epi; cells <- c(colnames(epi), paste0("t", seq_len(nrow(d$tme))))
  cnt <- Seurat::GetAssayData(epi, slot = "counts")
  meta <- data.frame(sample_id = c(epi$sample_id, d$tme$sample_id), tissue = c(epi$tissue, rep("NT", nrow(d$tme))),
                     major = rep(c("Epithelial", "Other"), c(ncol(epi), nrow(d$tme))),
                     cell_state = c(rep(NA_character_, ncol(epi)), d$tme$cell_state), row.names = cells)
  colnames(cnt) <- cells[seq_len(ncol(epi))]
  r <- suppressMessages(suppressWarnings(validate_hixscore(cnt, meta, verbose = FALSE)))
  expect_equal(r$layer, c("E_seHIX_z", "S_C_z", "S_M_z", "sysHIX"))
  expect_true(all(c("AUC_expected", "AUC", "pass") %in% names(r)))
  expect_false(isTRUE(attr(r, "all_pass")))                                                           # synthetic data is not the cohort
  expect_warning(suppressMessages(validate_hixscore(cnt, meta, verbose = FALSE)), "outside the expected values")
})
