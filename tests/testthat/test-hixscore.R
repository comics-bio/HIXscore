test_that("layer-by-layer runs with `previous` equal one ALL run", {
  skip_on_cran(); d <- synth_data()
  all_ <- suppressWarnings(suppressMessages(hixscore(d$epi, d$tme, layers = "ALL")))
  r <- suppressMessages(hixscore(d$epi, layers = "E"))
  expect_equal(r$status$status, c("done", "not run", "not run", "not run"))
  r <- suppressMessages(hixscore(d$epi, layers = "C", previous = r))
  expect_equal(r$status$status, c("reused", "done", "not run", "not run"))
  r <- suppressMessages(hixscore(d$epi, d$tme, layers = "M", previous = r))
  expect_equal(r$status$status, c("reused", "reused", "done", "not run"))
  r <- suppressMessages(hixscore(layers = "sysHIX", previous = r))             # no epi/tme needed: all layers reused
  expect_equal(r$status$status, c("reused", "reused", "reused", "done"))
  expect_equal(r$samples, all_$samples)
})

test_that("a failed layer does not stop the others and can be repaired", {
  skip_on_cran(); d <- synth_data()
  bad <- d$tme; bad$cell_state <- "not_a_state"                                 # M layer cannot be scored
  r <- suppressWarnings(suppressMessages(hixscore(d$epi, bad, layers = "ALL")))
  st <- setNames(r$status$status, r$status$layer)
  expect_equal(unname(st[c("E", "C", "M", "sysHIX")]), c("done", "done", "failed", "skipped"))
  expect_match(r$status$note[r$status$layer == "M"], "cell states")
  expect_match(r$status$note[r$status$layer == "sysHIX"], "missing layer")
  expect_true(all(r$samples$has_E) && all(r$samples$has_C) && !any(r$samples$has_M) && !any(r$samples$has_sysHIX))
  fixed <- suppressMessages(hixscore(d$epi, d$tme, layers = c("M", "sysHIX"), previous = r))
  expect_equal(unname(setNames(fixed$status$status, fixed$status$layer)[c("E", "C", "M", "sysHIX")]), c("reused", "reused", "done", "done"))
  expect_true(sum(fixed$samples$has_sysHIX) > 0)
  expect_error(suppressMessages(hixscore(d$epi, bad, layers = "ALL", on_error = "stop")), "Layer M failed")
})

test_that("layer arguments are validated and each layer asks for its own input", {
  skip_on_cran(); d <- synth_data()
  expect_error(hixscore(d$epi, layers = "X"), "Unknown layer")
  expect_equal(HIXscore:::.norm_layers(c("syshix", "e")), c("E", "sysHIX"))
  r <- suppressWarnings(suppressMessages(hixscore(d$epi, layers = "M")))          # M without tme
  expect_equal(r$status$status[r$status$layer == "M"], "failed"); expect_match(r$status$note[r$status$layer == "M"], "tme")
  r2 <- suppressWarnings(suppressMessages(hixscore(d$epi, layers = "sysHIX")))      # nothing to fuse
  expect_equal(r2$status$status[r2$status$layer == "sysHIX"], "skipped")
  expect_error(hixscore(d$epi, previous = list()), "previous")
})

test_that("M layer can use n_epi without epi", {
  skip_on_cran(); d <- synth_data(); n <- table(d$epi$sample_id)
  r <- suppressMessages(hixscore(tme = d$tme, layers = "M", n_epi = n))
  expect_equal(r$status$status[r$status$layer == "M"], "done"); expect_true(all(r$samples$has_M))
})

test_that("validate_MP_list checks names, gene vectors and size", {
  g <- paste0("G", 1:10)
  expect_error(validate_MP_list(list(A = g)), "at least 2")
  expect_error(validate_MP_list(list(g, g)), "names")
  expect_error(validate_MP_list(list(A = g, A = g)), "unique")
  expect_error(validate_MP_list(list(A = g, `1 bad` = g)), "syntactically valid")
  expect_error(validate_MP_list(list(A = g, B = 1:10)), "character")
  expect_error(validate_MP_list(list(A = g, B = g[1:3])), "fewer than")
  expect_warning(validate_MP_list(list(A = g, B = g), genes = g[1:3]), "Low gene coverage")
  expect_silent(validate_MP_list(default_MP_list()))
})

test_that("custom mp_list: mode 2 refits the E baseline on the user gene sets", {
  skip_on_cran(); d <- synth_data()
  mp <- default_MP_list()[1:4]; names(mp) <- paste0("Custom_", seq_along(mp))
  r <- suppressWarnings(suppressMessages(hixscore(d$epi, layers = "E", mode = "custom_baseline", mp_list = mp)))
  expect_equal(r$status$status[r$status$layer == "E"], "done")
  expect_true(all(r$samples$has_E))
  expect_gt(nrow(r$E), 0)
})

test_that("custom mp_list: mode 1 needs the shipped MP names", {
  skip_on_cran(); d <- synth_data()
  mp <- default_MP_list()[1:4]; names(mp) <- paste0("Custom_", seq_along(mp))
  expect_error(hixscore(d$epi, layers = "E", mode = "pretrained", mp_list = mp), "bound to the shipped MPs")
  r <- suppressWarnings(suppressMessages(hixscore(d$epi, layers = "E", mode = "pretrained", mp_list = default_MP_list())))
  expect_equal(r$status$status[r$status$layer == "E"], "done")
})

test_that("mp_list must contain the MPs of a supplied E_baseline", {
  skip_on_cran(); d <- synth_data()
  mp <- default_MP_list()[1:4]; names(mp) <- paste0("Custom_", seq_along(mp))
  tr <- suppressWarnings(suppressMessages(hixscore_train(d$epi, steps = "E", mp_list = mp)))
  expect_equal(tr$status$status[tr$status$layer == "E"], "done")
  expect_equal(names(tr$mp_list), names(mp))
  expect_equal(tr$E_baseline$mp_cols, names(mp))
  other <- default_MP_list()[5:8]; names(other) <- paste0("Other_", seq_along(other))
  expect_error(hixscore(d$epi, layers = "E", mode = "custom_baseline", E_baseline = tr$E_baseline, mp_list = other), "lacks MP")
  r <- suppressWarnings(suppressMessages(hixscore(d$epi, layers = "E", mode = "custom_baseline", E_baseline = tr$E_baseline, mp_list = mp)))
  expect_equal(r$status$status[r$status$layer == "E"], "done")
})

test_that("M_W builds the M-layer model from a user basis matrix", {
  skip_on_cran(); skip_if_not_installed("caret"); d <- synth_data()
  P <- build_M_composition(d$tme, "sample_id", "cell_state", min_cells = 200)
  W <- default_M_model()$W; W <- W[rownames(W) %in% rownames(P), , drop = FALSE]   # W states must be rows of the composition matrix
  tr <- suppressWarnings(suppressMessages(hixscore_train(tme = d$tme, steps = c("M_scan", "M"), M_W = W)))
  st <- setNames(tr$status$status, tr$status$layer)
  expect_equal(unname(st["M"]), "done"); expect_equal(unname(st["M_scan"]), "not run")
  expect_equal(rownames(tr$M_model$W), rownames(W)); expect_equal(ncol(tr$M_model$W), ncol(W))
  expect_error(hixscore_train(tme = d$tme, steps = "M", M_W = W, M_model = default_M_model()), "not both")
})
