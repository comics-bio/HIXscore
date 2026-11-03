# HIXscore 0.0.0.9000

* Three modes: `pretrained` (shipped oral model), `custom_baseline` (your NT samples), training (`hixscore_train()`).
* E, C and M layers are scored independently on the samples passing their own filters; samples are intersected only for sysHIX.
  The M layer keeps samples without epithelial cells and flags those whose epithelial cell number is below the training range.
* E-layer baseline fitting needs at least 3 NT samples (warning below 10).
* M-layer rank selection by the cophenetic rule rho(k-2) < rho(k-1) < rho(k) > rho(k+1).
* Cell annotation by two-level label transfer (`annotate_cells()`, `prepare_reference()`, `prepare_input()`); the published
  reference download is not available yet.
* Supports Seurat 4.x only.
