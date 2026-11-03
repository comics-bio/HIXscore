# Small synthetic data for tests: epithelial Seurat object + non-epithelial metadata
synth_data <- function(ns = 12, per = 45, seed = 11) {
  skip_if_not_installed("Seurat")
  set.seed(seed)
  mp <- default_MP_list(); pr <- default_BDM_programs(); st <- rownames(default_M_model()$W)
  genes <- unique(c(unlist(mp), unlist(pr), paste0("G", 1:2500))); n <- ns * per
  sid <- paste0("s", rep(seq_len(ns), each = per))
  tis <- rep(c("NT", "PD", "PrCA"), each = ns / 3)[as.integer(sub("s", "", sid))]
  mu <- rexp(length(genes), 0.5)
  m <- matrix(rpois(length(genes) * n, rep(mu, n)), length(genes), n, dimnames = list(genes, paste0("c", seq_len(n))))
  d <- unlist(c(pr, mp[1:3])); m[d, tis != "NT"] <- m[d, tis != "NT"] + rpois(length(d) * sum(tis != "NT"), 1)
  o <- suppressWarnings(Seurat::NormalizeData(Seurat::CreateSeuratObject(m), verbose = FALSE))
  o$sample_id <- sid; o$tissue <- tis; o$series <- paste0("ser", 1 + (as.integer(sub("s", "", sid)) %% 4))
  w <- default_M_model()$W
  tme <- data.frame(sample_id = rep(paste0("s", seq_len(ns)), each = 250), cell_state = NA_character_)
  tme$cell_state <- unlist(lapply(seq_len(ns), function(i) sample(st, 250, TRUE, prob = w[, 1 + i %% 4] + 0.01)))
  lab <- data.frame(sample_id = paste0("s", seq_len(ns)), label = rep(c("NT", "Inflammatory", "Malignant"), each = ns / 3))
  list(epi = o, tme = tme, labels = lab)
}
