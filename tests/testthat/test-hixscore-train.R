test_that("hixscore_train() steps run one by one and equal one ALL run", {
  skip_on_cran(); skip_if_not_installed("caret"); skip_if_not_installed("ordinalNet"); d <- synth_data()
  a <- suppressWarnings(suppressMessages(hixscore_train(d$epi, d$tme, d$labels, steps = "ALL")))
  st <- setNames(a$status$status, a$status$layer)
  expect_equal(unname(st[c("E", "C", "fusion")]), c("done", "done", "done"))
  expect_equal(unname(st[c("M_scan", "M", "BDM_scan", "BDM")]), rep("not run", 4))      # defaults: nothing derived
  b <- suppressWarnings(suppressMessages(hixscore_train(d$epi, steps = "E")))
  b <- suppressWarnings(suppressMessages(hixscore_train(d$epi, steps = "C", previous = b)))
  b <- suppressWarnings(suppressMessages(hixscore_train(d$epi, d$tme, d$labels, steps = "fusion", previous = b)))
  expect_equal(unname(setNames(b$status$status, b$status$layer)[c("E", "C", "fusion")]), c("reused", "reused", "done"))
  expect_equal(a$model$final_layer$beta, b$model$final_layer$beta)
  expect_equal(a$model$C_layer$beta, b$model$C_layer$beta)
})

test_that("the fusion step is skipped, not crashed, when a component is missing", {
  skip_on_cran(); d <- synth_data()
  r <- suppressWarnings(suppressMessages(hixscore_train(d$epi, d$tme, d$labels, steps = "fusion")))
  s <- r$status[r$status$layer == "fusion", ]
  expect_equal(s$status, "skipped"); expect_match(s$note, "missing component")
  expect_null(r$model)
  expect_error(hixscore_train(d$epi, d$tme, steps = "fusion"), "needs `labels`")
  expect_error(hixscore_train(d$epi, steps = "Z"), "Unknown step")
})

test_that("a failed step does not stop independent steps", {
  skip_on_cran(); d <- synth_data()
  e2 <- d$epi; e2$tissue <- "Tumor"                                              # no NT cells: E and C cannot be fitted
  r <- suppressWarnings(suppressMessages(hixscore_train(e2, d$tme, d$labels, steps = c("E", "C", "fusion"))))
  st <- setNames(r$status$status, r$status$layer)
  expect_equal(unname(st[c("E", "C")]), c("failed", "failed")); expect_equal(unname(st["fusion"]), "skipped")
  expect_match(r$status$note[r$status$layer == "fusion"], "failed/skipped")
  expect_error(suppressMessages(hixscore_train(e2, d$tme, d$labels, steps = "E", on_error = "stop")), "Step E failed|Layer E failed")
})

test_that("M steps derive the ecotypes and feed the fusion model", {
  skip_on_cran(); skip_if_not_installed("NMF"); skip_if_not_installed("caret"); d <- synth_data(ns = 12, per = 45)
  st <- rownames(default_M_model()$W)[1:20]; set.seed(5)
  tme <- data.frame(sample_id = rep(paste0("s", 1:12), each = 250), cell_state = unlist(lapply(1:12, function(i) sample(st, 250, TRUE, prob = c(rep(5, 5), rep(1, 15))[(seq_along(st) + i) %% 20 + 1]))))
  sc <- suppressWarnings(suppressMessages(hixscore_train(tme = tme, steps = "M_scan", M_model = "derive", M_ranks = 2:6)))
  expect_equal(sc$status$status[sc$status$layer == "M_scan"], "done"); expect_true(is.data.frame(sc$scans$M)); expect_true("candidate" %in% names(sc$scans$M))
  m <- suppressWarnings(suppressMessages(hixscore_train(tme = tme, steps = "M", M_model = "derive", M_k = 3, previous = sc)))
  expect_equal(ncol(m$M_model$W), 3); expect_equal(m$status$status[m$status$layer == "M_scan"], "reused")
  r <- suppressWarnings(suppressMessages(hixscore_train(d$epi, tme, d$labels, steps = c("E", "C", "fusion"), M_model = m$M_model, previous = m)))
  expect_equal(r$status$status[r$status$layer == "fusion"], "done"); expect_equal(ncol(r$M_model$W), 3)
  expect_true(any(grepl("M_CM3", r$model$M_layer$features_raw)))
})

test_that("BDM steps: bdm_k is needed for BDM, ALL stops after the scan", {
  skip_on_cran(); d <- synth_data()
  expect_error(hixscore_train(d$epi, steps = "BDM", programs = "derive"), "bdm_k")
  expect_error(hixscore_train(d$epi, steps = "BDM", programs = "derive", bdm_k = 3), NA)   # bdm_k given: reaches the BDM step (it is skipped without a scan)
  r <- suppressWarnings(suppressMessages(hixscore_train(d$epi, steps = "BDM", programs = "derive", bdm_k = 3)))
  expect_equal(r$status$status[r$status$layer == "BDM"], "skipped"); expect_match(r$status$note[r$status$layer == "BDM"], "BDM_scan")
  r2 <- suppressWarnings(hixscore_train(d$epi, d$tme, d$labels, steps = "ALL", programs = "derive", bdm_state_col = "tissue"))
  expect_equal(r2$status$status[r2$status$layer == "fusion"], "not run")                  # stopped before C/fusion, waiting for bdm_k
})

test_that("user supplied gene programs and M model are used as given", {
  skip_on_cran(); skip_if_not_installed("caret"); d <- synth_data()
  mine <- default_BDM_programs()[1:3]
  a <- suppressWarnings(suppressMessages(hixscore_train(d$epi, d$tme, d$labels, steps = c("E", "C", "fusion"), programs = mine)))
  expect_equal(names(a$programs), names(mine)); expect_equal(length(a$model$C_layer$norm_params), 4)
  expect_equal(length(a$C_baseline$nt_baseline), 3)
  expect_error(hixscore_train(d$epi, d$tme, d$labels, programs = list(A = 1:3, B = "x")), "character")
})

test_that("ALL with `previous` only runs the steps that are not finished", {
  skip_on_cran(); skip_if_not_installed("caret"); d <- synth_data()
  a <- suppressWarnings(suppressMessages(hixscore_train(d$epi, steps = c("E", "C"))))
  calls <- 0; real <- HIXscore:::fit_E_baseline
  testthat::local_mocked_bindings(fit_E_baseline = function(...) { calls <<- calls + 1; real(...) }, .package = "HIXscore")
  b <- suppressWarnings(suppressMessages(hixscore_train(d$epi, d$tme, d$labels, steps = "ALL", previous = a)))
  expect_equal(calls, 0)                                                         # E was finished: not fitted again
  expect_equal(unname(setNames(b$status$status, b$status$layer)[c("E", "C", "fusion")]), c("reused", "reused", "done"))
  suppressWarnings(suppressMessages(hixscore_train(d$epi, steps = "E", previous = a))); expect_equal(calls, 1)   # named step: rerun
})

test_that("bdm_args reach the BDM functions, including the pseudobulk settings", {
  skip_on_cran()
  calls <- list()
  testthat::local_mocked_bindings(
    bdm_constrained_genes = function(epi, ...) { calls$cg <<- list(...); list(genes = "g", frequency = NULL) },
    bdm_deviation_genes = function(epi, constrained, ...) { calls$dg <<- list(...); list(genes = "g", direction = c(g = "UP")) },
    bdm_full_deviation = function(epi, genes, ...) { calls$fd <<- list(...); matrix(0, 3, 3) },
    bdm_scan_modules = function(dev, ...) list(consensus = NULL, stability = data.frame(k = 2, mean_within_consensus = 1)),
    .package = "HIXscore")
  d <- synth_data(); d$epi$MP_state <- "MP_A"
  args <- list(max_var = 0.3, min_cells = 7, min_logcpm = 0.5, min_disease_types = 2, min_fc = 1.5, max_na = 0.4)
  suppressWarnings(suppressMessages(hixscore_train(d$epi, steps = "BDM_scan", programs = "derive", bdm_args = args)))
  expect_equal(calls$cg$max_var, 0.3); expect_equal(calls$cg$min_cells, 7); expect_equal(calls$cg$min_logcpm, 0.5)   # pseudobulk settings arrive
  expect_equal(calls$dg$min_disease_types, 2); expect_equal(calls$dg$min_fc, 1.5); expect_equal(calls$dg$min_cells, 7)
  expect_equal(calls$fd$max_na, 0.4); expect_equal(calls$fd$min_logcpm, 0.5)
  r <- suppressWarnings(suppressMessages(hixscore_train(d$epi, steps = "BDM_scan", programs = "derive", bdm_args = list(typo_arg = 1))))
  expect_equal(r$status$status[r$status$layer == "BDM_scan"], "failed"); expect_match(r$status$note[r$status$layer == "BDM_scan"], "Unknown `bdm_args`")
})

test_that("the bdm_args lists of hixscore_train() match the signatures of the BDM functions", {
  allowed <- function(f) setdiff(names(formals(f)), c("epi", "constrained", "genes", "dev", "...", "sample_col", "state_col", "tissue_col", "series_col", "nt_label"))
  src <- paste(deparse(HIXscore::hixscore_train), collapse = "\n")
  for (f in list(HIXscore::bdm_constrained_genes, HIXscore::bdm_deviation_genes, HIXscore::bdm_full_deviation, HIXscore::bdm_pseudobulk))
    for (a in allowed(f)) expect_true(grepl(paste0("\"", a, "\""), src, fixed = TRUE), info = paste("argument not allowed in bdm_args:", a))
})

test_that("leave-one-cohort-out needs at least 2 cohorts with enough samples", {
  meta <- data.frame(sample_id = paste0("s", 1:40), series = rep(c("a", "b", "c", "d"), each = 10))
  expect_equal(sort(HIXscore:::.loco_sets(meta)), c("a", "b", "c", "d"))
  expect_warning(HIXscore:::.loco_sets(meta[meta$series %in% c("a", "b"), ]), "Only 2 cohorts")
  expect_error(HIXscore:::.loco_sets(meta[meta$series == "a", ]), "at least 2 cohorts")
  few <- data.frame(sample_id = paste0("s", 1:12), series = rep(c("a", "b", "c"), each = 4))      # 4 samples each: none qualifies
  expect_error(HIXscore:::.loco_sets(few), "found 0")
})
