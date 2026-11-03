test_that("hixscore() pretrained: layers are scored independently, only sysHIX intersects", {
  skip_on_cran(); d <- synth_data()
  d$tme <- d$tme[!d$tme$sample_id %in% "s12", ]                                     # s12: no M layer
  d$tme <- rbind(d$tme, data.frame(sample_id = "extra", cell_state = sample(rownames(default_M_model()$W), 300, TRUE)))  # M only
  r <- suppressMessages(hixscore(d$epi, d$tme, mode = "pretrained"))
  s <- r$samples
  expect_true(all(c("E_seHIX", "E_seHIX_z", "S_C", "S_C_z", "S_M", "S_M_z", "sysHIX", "has_E", "has_C", "has_M", "has_sysHIX",
                    "n_epi_cells", "n_tme_cells", "M_epi_below_training_min") %in% names(s)))
  e <- s[s$sample_id == "extra", ]; expect_true(e$has_M && !e$has_E && !e$has_C && !e$has_sysHIX && is.na(e$sysHIX) && !is.na(e$S_M_z))
  expect_true(e$M_epi_below_training_min)                                          # no epithelial cells -> flagged
  x <- s[s$sample_id == "s12", ]; expect_true(x$has_E && x$has_C && !x$has_M && is.na(x$sysHIX) && !is.na(x$E_seHIX_z))
  expect_equal(sort(r$sysHIX$sample_id), sort(s$sample_id[s$has_sysHIX]))
  expect_false("s12" %in% r$sysHIX$sample_id)
  expect_equal(r$info$n_M_epi_below_training_min, 1)
})

test_that("hixscore() custom_baseline refits on NT and reports the NT sample count", {
  skip_on_cran(); d <- synth_data()
  r <- suppressWarnings(suppressMessages(hixscore(d$epi, d$tme, mode = "custom_baseline")))
  expect_equal(r$info$n_nt_E_baseline, 4)
  expect_true(sum(r$samples$has_sysHIX) > 0)
})

test_that("hixscore() custom_baseline with too few NT samples: E layer fails, the error is recorded", {
  skip_on_cran(); d <- synth_data()
  d$epi$tissue[d$epi$sample_id != "s1"] <- "PrCA"
  r <- suppressWarnings(suppressMessages(hixscore(d$epi, d$tme, mode = "custom_baseline")))
  e <- r$status[r$status$layer == "E", ]
  expect_equal(e$status, "failed"); expect_match(e$note, "NT sample")
  expect_equal(r$status$status[r$status$layer == "sysHIX"], "skipped")
  expect_error(suppressWarnings(suppressMessages(hixscore(d$epi, d$tme, mode = "custom_baseline", on_error = "stop"))), "NT sample")
})

test_that("custom programs need a model trained on them", {
  skip_on_cran(); d <- synth_data()
  expect_error(hixscore(d$epi, d$tme, mode = "custom_baseline", programs = default_BDM_programs()), "hixscore_train")
})

test_that("hixscore_train() works with three and two classes", {
  skip_on_cran(); skip_if_not_installed("caret"); skip_if_not_installed("ordinalNet"); d <- synth_data()
  for (lv in list(c("NT", "Inflammatory", "Malignant"), c("NT", "Disease"))) {
    lab <- d$labels; if (length(lv) == 2) lab$label <- ifelse(lab$label == "NT", "NT", "Disease")
    tr <- suppressWarnings(suppressMessages(hixscore_train(d$epi, d$tme, lab, levels = lv)))
    expect_s3_class(tr$model, "sysHIX_model")
    expect_equal(length(tr$model$metadata$severity_levels), length(lv))
    r <- suppressWarnings(suppressMessages(hixscore(d$epi, d$tme, mode = "custom_baseline", model = tr$model,
            E_baseline = tr$E_baseline, C_baseline = tr$C_baseline, M_model = tr$M_model)))
    v <- r$samples$sysHIX[r$samples$has_sysHIX]; expect_true(length(v) > 0)
    expect_true(all(v >= 0 & v <= length(lv) - 1))
  }
})

test_that("input errors are informative", {
  skip_on_cran(); d <- synth_data()
  d2 <- d$epi; d2$sample_id <- NULL
  expect_error(hixscore(d2, d$tme), "sample column")
  expect_error(hixscore(d$epi, d$tme[, c("sample_id")], mode = "pretrained"), "cell_state|columns")
})
