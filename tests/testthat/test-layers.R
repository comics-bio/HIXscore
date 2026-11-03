test_that("rank selection follows rho(k-2)<rho(k-1)<rho(k)>rho(k+1)", {
  expect_equal(select_M_rank(c(.90, .93, .97, .95, .92), 2:6)$optimal, 4)
  expect_true(is.na(select_M_rank(c(.99, .97, .95, .93), 2:5)$optimal))
  expect_true(is.na(select_M_rank(rep(1, 5), 2:6)$optimal))
})

test_that("M-layer scores are compositional and independent of other samples", {
  m <- default_M_model(); set.seed(1); st <- rownames(m$W)
  md <- data.frame(sample_id = rep(paste0("s", 1:4), each = 300), cell_state = sample(st, 1200, TRUE))
  all4 <- compute_M_CM_scores(md, m); two <- compute_M_CM_scores(md[md$sample_id %in% c("s1", "s2"), ], m)
  cm <- paste0("M_CM", 1:4)
  expect_equal(unname(rowSums(all4[, cm])), rep(1, 4), tolerance = 1e-8)
  expect_equal(all4[1:2, cm], two[, cm], ignore_attr = TRUE)
})

test_that("E-layer score with a frozen baseline does not depend on the other samples", {
  b <- readRDS(system.file("extdata", "E_seHIX_baseline_oral.rds", package = "HIXscore"))
  set.seed(2); n <- 1200
  md <- data.frame(sample_id = rep(paste0("s", 1:6), each = 200), tissue = "Tumor")
  for (k in b$mp_cols) md[[k]] <- rnorm(n)
  all6 <- compute_E_seHIX(md, baseline = b); two <- compute_E_seHIX(md[md$sample_id %in% c("s1", "s2"), ], baseline = b)
  expect_equal(all6$E_seHIX[1:2], two$E_seHIX)
})

test_that("shipped final layer maps z-scores to a sysHIX in [0, 2]", {
  skip_if_not_installed("ordinalNet"); requireNamespace("ordinalNet")
  m <- default_sysHIX_model()
  X <- matrix(c(-2, 0, 3, -2, 0, 3, -2, 0, 3), 3, dimnames = list(NULL, m$final_layer$features))
  cum <- predict(m$final_layer$model, newx = X, type = "response", lambda = m$final_layer$lambda)
  s <- as.numeric(HIXscore:::.cum_to_class(cum, m$metadata$severity_levels) %*% m$metadata$severity_mapping)
  expect_true(all(s >= 0 & s <= 2)); expect_true(s[3] > s[1])
})

test_that("BDM programs are validated", {
  expect_error(validate_BDM_programs(list(A = 1:3, B = "x")))
  expect_error(validate_BDM_programs(list(A = letters, A = letters)))
  expect_silent(validate_BDM_programs(default_BDM_programs()))
})

test_that("evaluate_scores gives AUC 1 for perfectly separated groups", {
  s <- data.frame(tissue = rep(c("NT", "D"), each = 6), sysHIX = c(1:6, 11:16))
  ev <- evaluate_scores(s, scores = "sysHIX", stages = c("NT", "D"))
  expect_equal(ev$auc_nt_vs_disease$auc, 1)
})

test_that("layers are scored independently and only sysHIX intersects samples", {
  m <- default_sysHIX_model(); mm <- default_M_model(); set.seed(3)
  st <- rownames(mm$W)
  tme <- data.frame(sample_id = rep(c("a", "b", "c"), each = 300), cell_state = sample(st, 900, TRUE))
  M <- compute_M_CM_scores(tme, mm)
  n_epi <- c(a = 100, b = 0)                      # b has no epithelial cells, c is unknown
  Ms <- compute_S_M(M, n_epi, m)
  expect_equal(sort(Ms$sample_id), c("a", "b", "c"))                # nobody dropped
  expect_equal(Ms$M_epi_below_training_min[match(c("a", "b", "c"), Ms$sample_id)], c(FALSE, TRUE, TRUE))
  expect_false(anyNA(Ms$S_M_z))
  C <- data.frame(sample_id = c("a", "x"), n_epi_cells = c(100, 80), C_Low_Frac = .3, C_Mid_Frac = .3, C_High_Frac = .4, C_mean = .4)
  C <- HIXscore:::.S_C_from_features(C, m)
  E <- data.frame(sample_id = c("a", "x", "y"), n_cells = c(100, 80, 60), E_seHIX = c(.5, .6, .7))
  E <- compute_E_z(E, m)
  S <- predict_sysHIX(E, C, Ms, m)
  expect_equal(S$sample_id, "a")                                    # intersection only here
  expect_true(all(c("a", "b", "c") %in% Ms$sample_id) && "y" %in% E$sample_id)
})

test_that("E baseline needs enough NT samples", {
  mp <- paste0("MP_", 1:12, "_X"); set.seed(4)
  mk <- function(n_nt) { n <- (n_nt + 6) * 40
    md <- data.frame(sample_id = rep(paste0("s", seq_len(n_nt + 6)), each = 40))
    md$tissue <- ifelse(as.integer(sub("s", "", md$sample_id)) <= n_nt, "NT", "Tumor")
    for (k in mp) md[[k]] <- rnorm(n); md }
  expect_error(fit_E_baseline(mk(2), min_cells = 30), "Only 2 NT sample")
  expect_warning(b <- fit_E_baseline(mk(4), min_cells = 30), "unstable")
  expect_equal(b$n_nt, 4)
  expect_no_warning(fit_E_baseline(mk(10), min_cells = 30))
})

test_that("evaluate_scores can exclude samples with an extrapolated M-layer adjustment", {
  s <- data.frame(tissue = rep(c("NT", "D"), each = 8), S_M_z = c(1:8, 9:16),
                  M_epi_below_training_min = c(rep(FALSE, 14), TRUE, TRUE))
  a <- evaluate_scores(s, scores = "S_M_z", stages = c("NT", "D"))
  b <- evaluate_scores(s, scores = "S_M_z", stages = c("NT", "D"), exclude_flagged = TRUE)
  expect_equal(a$auc_nt_vs_disease$n_hi, 8); expect_equal(b$auc_nt_vs_disease$n_hi, 6)
})

test_that("Seurat >= 5 is refused with a clear message", {
  skip_if(utils::packageVersion("SeuratObject") >= "5.0.0")
  expect_true(HIXscore:::.check_seurat())
  testthat::local_mocked_bindings(packageVersion = function(pkg, ...) package_version("5.1.0"), .package = "utils")
  expect_error(HIXscore:::.check_seurat(), "supports Seurat 4.x only")
})

test_that(".need names the missing package", {
  expect_error(HIXscore:::.need("notARealPackageXYZ", "testing"), "notARealPackageXYZ.*required for testing")
})

test_that("M-layer reproduces the shipped training scores (regression)", {
  m <- default_M_model(); V <- m$V_train; H <- m$H_train
  H2 <- vapply(seq_len(ncol(V)), function(i) { h <- nnls::nnls(m$W, V[, i])$x; h / sum(h) }, numeric(nrow(H)))
  expect_equal(unname(H2), unname(H), tolerance = 1e-8)
})
