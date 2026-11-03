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
