# Rebuild inst/extdata from the training reference parameters (build_sysHIX_reference_params.R output)
root <- "/data/comics-zhangyf/2025.06_oral_homeostasis"
ref  <- readRDS(file.path(root, "step_5_sysHIX/sysHIX_reference_params.rds"))
ext  <- "inst/extdata"

# MP gene lists (current MP naming used with the reference parameters)
e <- new.env(); load(file.path(root, "step_2_epi_NMF/MP_list_named_new.rdata"), e)
mp <- e$MP_list_named
stopifnot(all(ref$E$mp_cols %in% names(mp)))
saveRDS(mp, file.path(ext, "MP_list.rds"))

# E layer: frozen per-feature mean/SD and NT baseline of the training set
E <- ref$E
saveRDS(structure(list(mp_cols = E$mp_cols, features = E$feature_cols,
                       mean = E$feat_mean, sd = E$feat_sd, nt_baseline = unname(E$nt_baseline),
                       n_nt = E$n_nt_samples, n_train_samples = E$n_train_samples,
                       min_cells = 10, creation_date = as.Date(ref$created)),
                  class = "sysHIX_E_baseline"), file.path(ext, "E_seHIX_baseline_oral.rds"))

# C layer: BDM programs and training NT baseline
saveRDS(ref$scHIX$baseline_programs, file.path(ext, "BDM_programs.rds"))
saveRDS(structure(list(nt_baseline = ref$scHIX$nt_baseline, programs = ref$scHIX$baseline_programs,
                       n_nt_cells = NA_integer_, creation_date = as.Date(ref$created)),
                  class = "sysHIX_C_baseline"), file.path(ext, "scHIX_baseline_oral.rds"))

# M layer: NMF basis and frozen training composition scaling
CM <- ref$CM; W <- CM$W; rownames(W) <- CM$celltypes; colnames(W) <- paste0("CM", seq_len(ncol(W)))
saveRDS(list(W = W, comp_min = CM$row_min[CM$celltypes],
             comp_max = (CM$row_min + CM$row_range)[CM$celltypes], k = ncol(W),
             V_train = CM$V_train, H_train = CM$H_nnls_train,
             min_total_cells = CM$min_total_cells, creation_date = as.Date(ref$created)),
        file.path(ext, "M_model_oral_K4.rds"))

file.copy(file.path(root, "step_5_sysHIX/sysHIX_model_complete_new.rds"),
          file.path(ext, "sysHIX_model_oral.rds"), overwrite = TRUE)
