# Build the public validation data set for HIXscore::validate_hixscore() from the annotated oral independent cohorts.
# The cohorts, label harmonisation and gene filter are those of sysHIX_run_Oral_public_new3.R.
#
# Output (written to `out_dir`):
#   validation_counts.rds  raw counts (dgCMatrix, genes x epithelial cells); cells are numbered 1..N, no barcodes
#   validation_meta.rds    one row per cell (epithelial cells first, same numbers as the count columns):
#                          sample_id (anonymous, S01...), tissue, major, cell_state
# No cohort, study, patient or barcode information is kept.
#
# Usage: Rscript make_validation_data.R [out_dir] [lib_path] [prep_rds]
#   prep_rds: optional cache with list(counts = <epithelial counts after the gene filter>, meta = <all cells>), written by
#   the first part of this script; when it exists the 13 GB cohort files are not loaded again.
args <- commandArgs(TRUE)
out_dir <- if (length(args) >= 1) args[1] else "."
if (length(args) >= 2 && nzchar(args[2])) .libPaths(c(args[2], .libPaths()))
prep_rds <- if (length(args) >= 3) args[3] else file.path(out_dir, "prep_oral_independent.rds")
suppressMessages({library(Seurat); library(Matrix)})

if (file.exists(prep_rds)) {
  prep <- readRDS(prep_rds)
} else {
  d <- "~/data/2025.06_oral_homeostasis/step_5_sysHIX/"
  for (f in c("oral_scRNA_hra340", "oral_scRNA_hra370", "oral_scRNA_prj396", "oral_scrna_gse789_oral", "oral_scRNA_gse630"))
    load(paste0(d, f, "_anno_transfer.rdata"))
  scRNA_hra340_anno$condition[scRNA_hra340_anno$condition == "Normal"] <- "NT"
  scRNA_hra370_anno$condition[scRNA_hra370_anno$condition == "Healthy"] <- "NT"
  scRNA_prj396_anno$condition[scRNA_prj396_anno$condition == "HC"] <- "NT"
  scrna_gse789_oral_anno$condition <- scrna_gse789_oral_anno$tissue
  scrna_gse789_oral_anno$condition[scrna_gse789_oral_anno$condition == "LP"] <- "OLP"
  scRNA_gse630_anno$condition <- scRNA_gse630_anno$tissue
  scRNA_hra340_anno$sample_id <- scRNA_hra340_anno$sample_origin
  scRNA_hra370_anno$sample_id <- scRNA_hra370_anno$sample
  nms <- c("scRNA_hra340_anno", "scRNA_hra370_anno", "scRNA_prj396_anno", "scrna_gse789_oral_anno", "scRNA_gse630_anno")
  sel <- Reduce(intersect, lapply(nms, function(n) colnames(get(n)@meta.data)))
  for (n in nms) { o <- get(n); o@meta.data <- o@meta.data[, sel]; assign(n, o) }
  all_obj <- merge(scRNA_hra340_anno, list(scRNA_hra370_anno, scRNA_prj396_anno, scrna_gse789_oral_anno, scRNA_gse630_anno))
  rm(list = nms); gc()
  all_obj$tissue <- all_obj$condition
  cnt <- GetAssayData(all_obj, slot = "counts")
  all_obj <- subset(all_obj, features = rownames(cnt)[rowSums(cnt > 0) >= 30]); rm(cnt)
  is_epi <- all_obj$predicted.major == "Epithelial"
  prep <- list(counts = GetAssayData(all_obj[, is_epi], slot = "counts"), meta = all_obj@meta.data)
  saveRDS(prep, prep_rds)
}

meta <- prep$meta; counts <- prep$counts
is_epi <- meta$predicted.major == "Epithelial"
stopifnot(identical(colnames(counts), rownames(meta)[is_epi]), all(meta$tissue %in% c("NT", "OSF", "CHC", "OLK", "OLP", "Tumor")))

# anonymous sample ids in random order (no link to cohort or patient, not sorted by cohort)
set.seed(20261126)
sid <- unique(as.character(meta$sample_id))
map <- setNames(sprintf("S%02d", seq_along(sid)), sample(sid))
ord <- c(which(is_epi), which(!is_epi))                  # epithelial cells first, as the count columns
out_meta <- data.frame(sample_id = unname(map[as.character(meta$sample_id)])[ord],
                       tissue = as.character(meta$tissue)[ord],
                       major = as.character(meta$predicted.major)[ord],
                       cell_state = as.character(meta$predicted.minor)[ord],
                       stringsAsFactors = FALSE)
rownames(out_meta) <- as.character(seq_len(nrow(out_meta)))
colnames(counts) <- as.character(seq_len(ncol(counts)))
stopifnot(identical(colnames(counts), rownames(out_meta)[out_meta$major == "Epithelial"]))
counts <- as(counts, "CsparseMatrix")

saveRDS(counts, file.path(out_dir, "validation_counts.rds"), compress = "xz")
saveRDS(out_meta, file.path(out_dir, "validation_meta.rds"), compress = "xz")
for (f in c("validation_counts.rds", "validation_meta.rds")) {
  p <- file.path(out_dir, f)
  cat(sprintf("%-24s size=%d md5=%s\n", f, file.size(p), tools::md5sum(p)))
}
cat(sprintf("epithelial cells=%d, all cells=%d, genes=%d, samples=%d\n", ncol(counts), nrow(out_meta), nrow(counts), length(map)))
print(table(out_meta$tissue)); print(table(out_meta$major))
