# HIXscore

Multi-layer epithelial instability scores from single-cell RNA-seq, fused into **sysHIX**:

| Layer | Score | Input cells | What it measures |
|---|---|---|---|
| E | seHIX | epithelial | distance of the sample's meta-program profile (activity + heterogeneity) to the normal-tissue (NT) profile |
| C | scHIX -> S_C | epithelial | cell-level distance of baseline-deviation-module (BDM) scores to the NT baseline, summarised per sample |
| M | CM -> S_M | immune + stromal | cell-state composition projected onto NMF ecotypes (CMs) |
| fusion | sysHIX | all | ordinal-logistic fusion of E, S_C and S_M; expected class index (0-2 with three classes) |

All scores are **sample-level**. The package is a research tool, not a diagnostic.

## Install

```r
# install.packages("remotes")
remotes::install_github("comics-bio/HIXscore")
```

**Supported Seurat versions:** Seurat 4.3 - 4.4. Seurat >= 5 is not supported; the main functions stop with an
explanatory error.

## Workflow

```
0. reference  ->  1. annotation  ->  2. hixscore()  [mode 1 / mode 2]
                                  ->  2. hixscore_train() -> hixscore()  [mode 3]
                                  ->  3. evaluate_scores() / plot_scores()
```

Every mode starts from cells with a **major cell type** and a **cell state**. The E and C layers use the epithelial
cells, the M layer uses the immune and stromal cells and needs the cell-state labels of the M-layer model (the 72 states
of `rownames(default_M_model()$W)`, e.g. `Mac_MP_1_M2_like`). If your cells carry no such labels, annotate them first
(steps 0 and 1).

## 0. Reference data

The annotation transfers labels from a single-cell reference.

### Use the HIXscore reference (default)

```r
library(HIXscore)
ref <- load_reference()                  # downloads ~550 MB once into the user cache, checks its md5 sum, prepares it
```

`load_reference()` downloads `scRNA_ref_HIXscore.rds` from
<https://comics.med.sustech.edu.cn/track/scRNA-ref/scRNA_ref_HIXscore.rds> into `tools::R_user_dir("HIXscore", "cache")`
and reuses the cached copy afterwards. To download it yourself, use `download_reference()` or any other tool and pass the
file: `load_reference(path = "scRNA_ref_HIXscore.rds")`.

Preparing the reference (splitting it by major cell type and recomputing a PCA within each type) is done once per R
session. Save the prepared reference to skip it next time:

```r
save_reference(ref, "HIXscore_reference_prepared.rds")
ref <- load_reference(path = "HIXscore_reference_prepared.rds")   # prepared references are loaded as they are
```

Reading an `.rds` file runs R code stored in the file, and the md5 sum only detects damaged or changed downloads. Load
files from sources you trust.

### Use your own reference

```r
ref <- load_reference(path = "my_reference.rds",
                      major_col = "celltype_major", minor_col = "cell_type2", reduction = "pca",
                      ref_features = NULL)       # optional marker genes per major type, see below
```

A reference of your own has to contain:

* a **Seurat (v4) object** with the **raw counts** in the `RNA` assay (all cells, gene symbols as row names);
* a metadata column with the **major cell type** (`major_col`) and one with the **cell state** (`minor_col`), without
  missing values. Epithelial cells must be labelled `"Epithelial"` in the major-type column, and every cell state belongs
  to one major type;
* a **PCA reduction** (`reduction`) with the embeddings of **all** cells and the **feature loadings** (the genes used
  for the PCA, usually the variable genes);
* at least 50 cells per cell state and at least 100 cells per major type (smaller ones are not used for the second level).

`ref_features` is an optional named list (names = major types) of marker genes of their cell states; they are added to
the genes used to transfer the cell states. `default_ref_features()` is the list belonging to the HIXscore reference.

For **modes 1 and 2** the cell states of the reference must include the 72 immune and stromal states of the shipped
M-layer model; `check_reference(ref)` lists the missing ones (`load_reference()` warns about them). A reference with a
different set of cell states can only be used together with your own M layer (mode 3).

## 1. Annotate the cells

```r
# obj: Seurat object with raw counts, all cells of all samples, unique cell names, a sample column
obj <- annotate_cells(obj, ref, sample_col = "sample_id",
                      checkpoint_dir = "annotation_checkpoints")
attr(obj, "annotation_summary")         # annotated fraction, cells per major type, fraction in the model's states
attr(obj, "annotation_failed")          # samples that could not be annotated, with the reason
```

* Each sample is annotated on its own, so the labels of a sample do not depend on the other samples.
* A sample that fails does not stop the others. With `checkpoint_dir` every finished sample is stored, and a repeated
  call annotates only the missing samples. Use a new directory when you change the reference or the settings.
* Samples with fewer than `min_cells_per_sample` (50) cells are not annotated. A major type with fewer than
  `min_query_cells` (50) cells in a sample gets no cell state (`NA`).
* By default the cell state is predicted for the non-epithelial cells only (`minor_for = "non-epithelial"`). Use
  `minor_for = "all"` when you need the epithelial states, which the BDM steps of mode 3 use.
* The reduction of the reference is used as prepared (`pca` by default); `reduction = ` selects another one when the
  reference was prepared with it.

Then split the cells into the two inputs of `hixscore()`:

```r
inp <- prepare_input(obj)               # genes detected in < 30 cells are dropped, cells are log-normalised together
epi <- inp$epi                          # epithelial cells (Seurat)   -> E and C layers
tme <- inp$tme                          # non-epithelial cells        -> M layer (sample_id, cell_state)
inp$summary
```

`prepare_input(obj, annotate = TRUE, reference = ref)` runs the annotation as well. If your cells already carry
`predicted.major` / `predicted.minor`, or other column names (`major_col`, `minor_col`), `prepare_input()` uses them
directly; the cell states have to follow the vocabulary of the M-layer model.

Input requirements: `epi` and `tme` have a sample column (`sample_id`); for mode 2 `epi` also needs the tissue column
(`tissue`, with the NT label `"NT"`), for mode 3 also a `series` column (study/batch). At least 30 epithelial cells and
200 immune/stromal cells per sample are used by default (`min_cells_epi`, `min_cells_tme`).

## 2. Run the layers: `hixscore()` (modes 1 and 2)

```r
res <- hixscore(epi, tme, layers = "ALL", mode = "pretrained")
res$samples        # one row per sample: E, S_C, S_M, their z-scores, sysHIX, has_E/has_C/has_M/has_sysHIX
res$status         # state, time and error message of every layer
```

`layers` selects what is computed: any of `"E"`, `"C"`, `"M"`, `"sysHIX"`, or `"ALL"`. Each layer needs only its own input
and is scored on the samples that pass its own filters; samples are intersected only in the `sysHIX` step. A layer that
fails does not stop the others: the error is stored in `res$status`, and `sysHIX` is skipped when a layer is missing.
Pass the result as `previous` to run only what is missing:

```r
res <- hixscore(epi, tme, layers = "E")                          # run the layers one at a time ...
res <- hixscore(epi, tme, layers = "C", previous = res)
res <- hixscore(epi, tme, layers = "M", previous = res)          # ... M needs tme and the epithelial cell numbers (epi or n_epi)
res <- hixscore(layers = "sysHIX", previous = res)               # fusion of the finished layers

res <- hixscore(epi, tme, layers = c("M", "sysHIX"), previous = res)   # repeat one layer and the fusion
```

Layers that are not requested are taken from `previous` (status `reused`). Splitting the run does not make it shorter; it
lets you control the steps and repeat only the one that failed. The time of each layer is in `res$status`.

* **Mode 1** (`mode = "pretrained"`): shipped oral baselines and model.
* **Mode 2** (`mode = "custom_baseline"`): the E and C baselines and the scHIX tertile thresholds are refitted on the NT
  samples of `epi`; the layer-to-score mappings and the fusion model stay at the oral training values.

```r
res <- hixscore(epi, tme, mode = "custom_baseline", tissue_col = "tissue", nt_label = "NT")
```

**NT sample number (mode 2).** The E-layer baseline and feature standardisation are estimated from your NT samples:
fewer than 3 NT samples (with enough epithelial cells) is an error and 3 to 9 gives a warning. The number used is stored in
`res$info$n_nt_E_baseline`.

**Samples with few epithelial cells.** The M layer keeps samples that have no or few epithelial cells. The model removes
the effect of the epithelial cell number from every CM with a regression fitted on training samples with more than 10
epithelial cells; samples with <= 10 epithelial cells are flagged (`M_epi_below_training_min`, counted in `res$info`).
`min_epi_M` sets a minimum number of epithelial cells for the M layer, and
`evaluate_scores(..., exclude_flagged = TRUE)` leaves the flagged samples out.

## 3. Train your own model: `hixscore_train()` (mode 3)

`hixscore_train()` runs the training steps one at a time or together:

| `steps` | What it does | Result |
|---|---|---|
| `"M_scan"` | NMF rank scan of the cell-state composition | `$scans$M` (cophenetic correlation per rank, `candidate`, `optimal_rank`) |
| `"M"` | cell-state ecotypes (CMs) at rank `M_k` | `$M_model` |
| `"BDM_scan"` | leave-one-cohort-out constrained and deviating genes, consensus-clustering scan | `$scans$bdm$stability`, `$intermediate` |
| `"BDM"` | BDM gene programs at `bdm_k` | `$programs` |
| `"E"` | E-layer baseline from the NT samples | `$E_baseline` |
| `"C"` | C-layer baseline and scHIX thresholds from the NT cells | `$C_baseline` |
| `"fusion"` | layer scores, out-of-fold S_C and S_M, fusion model | `$model` |
| `"ALL"` | every step that is needed | |

`M_model` and `programs` decide where the M-layer ecotypes and the BDM gene programs come from: `"default"` (shipped),
`"derive"` (computed from your data in the `M*` / `BDM*` steps) or an object you supply (a model from
`build_M_ecotypes()`, or a named list of gene sets). Only `"derive"` runs the corresponding steps.

```r
labels <- data.frame(sample_id = ..., label = ...)    # label: NT / Inflammatory / Malignant, or NT / Disease

# (a) only the E, C and fusion steps, shipped M layer and BDM programs
tr <- hixscore_train(epi, tme, labels, steps = c("E", "C", "fusion"),
                     levels = c("NT", "Inflammatory", "Malignant"))   # levels = c("NT", "Disease") for two classes

# (b) everything, with the steps run one by one and finished steps reused
tr <- hixscore_train(epi, tme, labels, steps = "M_scan", M_model = "derive")
tr$scans$M                                             # inspect, then choose M_k (or leave "auto")
tr <- hixscore_train(epi, tme, labels, steps = "M", M_model = "derive", M_k = 4, previous = tr)
```

Any step can fail without stopping the others; `tr$status` lists the state, time and error of each step, and `previous`
repeats only what is missing. With `previous`, `steps = "ALL"` runs the unfinished steps only; a step named explicitly is
always run again.

Use the result with `hixscore()`:

```r
res <- hixscore(new_epi, new_tme, mode = "custom_baseline", model = tr$model, E_baseline = tr$E_baseline,
                C_baseline = tr$C_baseline, M_model = tr$M_model, programs = tr$programs)
```

### M layer: ecotypes from your data

`M_model = "derive"` builds the composition matrix of your cell states (`build_M_composition()`), scans the NMF rank
(`scan_M_rank()`) and builds the ecotypes (`build_M_ecotypes()`). `M_k = "auto"` takes the rank that satisfies
rho(k-2) < rho(k-1) < rho(k) and rho(k) > rho(k+1) for the cophenetic correlation rho; the scan has to cover ranks on both
sides, and no rank is returned if none satisfies the rule. Your own NMF basis `W` (states x CMs, row names = states) can be
used with `build_M_ecotypes(P, W = W)`, and the resulting model passed as `M_model`.

### C layer: selecting the BDM gene sets

The BDM gene programs are selected from your data in two steps because the number of modules has to be chosen from the
stability table.

**Input.** The epithelial cells `epi` need the raw counts (`counts` slot), the sample column, a tissue column with NT
samples and disease samples, a `series` column (study/batch), and the **epithelial state of every cell**
(`bdm_state_col`). The selection is repeated leaving each cohort out: it uses the cohorts with more than 5 samples and
stops with an error when fewer than 2 qualify (2 cohorts give a warning because the frequency across iterations is then
coarse). A gene has to deviate in at least `min_disease_types` disease types (default 3; when your data have fewer disease
types, all of them are required). The epithelial states come from the annotation with `minor_for = "all"`; use
`bdm_state_col = "predicted.minor"` for the labels of `annotate_cells()`.

```r
# annotate with the epithelial states
obj <- annotate_cells(obj, ref, minor_for = "all", checkpoint_dir = "annotation_checkpoints")

# step 1: constrained genes, deviating genes and the consensus-clustering scan (the slow part)
tr <- hixscore_train(epi, steps = "BDM_scan", programs = "derive", bdm_state_col = "predicted.minor",
                     bdm_args = list(), bdm_maxK = 10, bdm_reps = 200)
tr$scans$bdm$stability          # mean within-module consensus for k = 2 ... bdm_maxK
tr$intermediate$deviation$genes # the genes entering the clustering (direction in $direction)

# step 2: gene programs at the chosen number of modules (takes seconds, the scan is reused)
tr <- hixscore_train(epi, steps = "BDM", programs = "derive", bdm_k = 5, bdm_drop = integer(), previous = tr)
tr$programs                     # named list Module_1 ... of gene symbols
```

What the steps do: (1) **constrained genes** - per leave-one-cohort-out iteration, `variancePartition` on the NT
pseudobulk (sample x epithelial state) keeps genes for which state, sample, subsite and cohort each explain less than
`max_var` (0.15) of the variance; genes kept in at least `min_freq` (0.7) of the iterations are retained. (2)
**deviating genes** - per iteration and disease type, genes whose deviation from the NT mean of the same state exceeds
`min_fc` (1.2-fold) in the same direction in at least `min_consistency` (0.7) of the pseudobulk samples; genes needed in at
least `min_disease_types` disease types and in `min_freq` of the iterations. (3) **module scan** - consensus clustering of
the genes by their co-deviation for k = 2 to `bdm_maxK`. (4) **programs** - for the chosen k, each module keeps the genes
of its majority direction (up or down); `bdm_drop` lists cluster ids to discard after inspection.

Arguments of the first three steps (`max_var`, `min_freq`, `min_fc`, `min_consistency`, `min_disease_types`,
`min_cells`, `min_logcpm`, ...) are given in `bdm_args`; see `?bdm_constrained_genes`, `?bdm_deviation_genes`,
`?bdm_pseudobulk`. The steps can also be called directly:

```r
cg  <- bdm_constrained_genes(epi, state_col = "predicted.minor")
dg  <- bdm_deviation_genes(epi, cg$genes, state_col = "predicted.minor")
dev <- bdm_full_deviation(epi, dg$genes, state_col = "predicted.minor")
sc  <- bdm_scan_modules(dev, maxK = 10)
programs <- build_BDM_geneset(dev, dg$direction, k = 5, scan = sc)
```

`steps = "ALL"` with `programs = "derive"` and no `bdm_k` stops after the scan and tells you to choose `bdm_k`. **Your own
gene sets** are used with `programs = list(Module_1 = c(...), Module_2 = c(...), ...)` (at least two sets with at least
five genes each, checked by `validate_BDM_programs()`); a model trained with them has to be used with the same gene sets.

### Notes on the trained model

* S_C and S_M of the training samples are **out-of-fold** (folds grouped by `series`); the fusion model itself is fitted on
  those, and the returned training sysHIX is computed on the samples the fusion model was fitted on.
* Two classes: the fusion is a penalised logistic regression and sysHIX is the **probability of disease** (0-1), not the 0-2
  expected class index of three classes.
* A model trained with custom gene sets cannot be combined with the shipped model.

## 4. Evaluate

```r
ev <- evaluate_scores(res$samples, stages = c(NT = "NT", OLP = "non-malignant", Tumor = "malignant"))
plot_scores(res$samples, stages = c(NT = "NT", OLP = "non-malignant", Tumor = "malignant"))
```

`evaluate_scores()` gives the AUC of NT against the other groups (DeLong CI), pairwise comparisons (Wilcoxon, AUC, BH
adjustment) and, for three or more ordered groups, Spearman and linear-trend tests. The AUC direction is fixed so that a
higher stage has a higher score. Groups with <= 3 samples are skipped.

## Notes

* `Seurat::AddModuleScore` selects control genes from the cells it is given, so the module scores of a sample depend on
  which cells are scored together. Samples below the cell threshold are removed before the module scoring.
* In mode 2 only the E and C baselines and the scHIX thresholds are refitted; the downstream mappings are not.
