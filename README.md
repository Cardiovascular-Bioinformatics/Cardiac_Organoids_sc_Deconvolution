# RCTD deconvolution of vascularised cardiac organoid (VCO) bulk RNA-seq

Analysis code for the cell-type deconvolution reported in Fig. 4G–J (reference atlas: Supplementary Fig. 5A–C) of:

> Kuburas R, *et al.* TINF2 Drives Mitophagy Impairment and Mitochondrial Dysfunction in Type 2 Diabetes Revealed by Donor-Derived Vascularised Cardiac Organoids. [Journal, year, DOI]

## Overview

Bulk RNA-seq profiles of VCOs from three donors with type 2 diabetes (DB VCO; RK001–RK003) and three non-diabetic donors (ND VCO; RK004–RK006) were deconvolved with RCTD (R package `spacexr`). The reference was a re-analysed and annotated cardiac organoid scRNA-seq atlas (GEO: GSE263193).

RCTD is run in **full mode** (`doublet_mode = "full"`). Each bulk library is treated as one pixel. The per-sample weights (`rctd@results$weights`) are converted to proportions with `spacexr::normalize_weights()`, so they sum to 1 in each sample. **Every panel of Fig. 4G–J uses this single matrix of normalised proportions** (`results/02_full_mode_normalised_proportions.csv`).

## Scripts

Run from the repository root, in order:

| Script | Purpose | Main outputs |
|---|---|---|
| `01_reference_atlas_GSE263193.R` | QC, clustering (Louvain, resolution 0.2) and annotation of GSE263193 | Annotated Seurat object; Supp. Fig. 5A–C |
| `02_rctd_full_mode.R` | RCTD full mode; normalisation to proportions | Raw weights, normalised proportions, QC, genes used |
| `03_original_vs_rerun.R` | (Optional) Reproduces the values plotted in the previous figure version and compares them with the rerun | Original weights, concordance |
| `04_figures.R` | Fig. 4G–J, per-group summaries and Welch t-tests | Fig. 4G–J; group summary; statistics |

```r
source("01_reference_atlas_GSE263193.R")
source("02_rctd_full_mode.R")
source("03_original_vs_rerun.R")   # optional; requires the original RCTD object
source("04_figures.R")
```

## Key settings

| Setting | Value |
|---|---|
| RCTD mode | `doublet_mode = "full"` |
| Reference | raw UMI counts, 12 annotated cell types, `n_max_cells = 10000`, `CELL_MIN_INSTANCE = 25` |
| Bulk input | raw read counts, duplicated gene symbols summed, `nUMI` = total counts per library |
| Filtering | `UMI_min = 100`, `UMI_max = 1e9` (so that no bulk library is excluded) |
| Gene selection | spacexr defaults (`gene_cutoff_reg = 2e-4`, `fc_cutoff_reg = 0.75`) |
| Output used | `rctd@results$weights`, normalised per sample with `normalize_weights()` |
| Random seed | `20261102` |
| Statistics | two-sided Welch's t-test, DB vs ND VCO, n = 3 donors per group |


## Inputs

Bulk RNA-seq will be made available upon request after the publication of the manuscript.

## Outputs included in this repository

Result files contains the RCTD outputs used for Fig. 4G–J, so they can be inspected without re-running the analysis:

| File | Content |
|---|---|
| `02_full_mode_raw_weights.csv` | `rctd@results$weights` (full mode, not constrained to sum to 1) |
| `02_full_mode_normalised_proportions.csv` | Raw weights divided by their per-sample sum; used for all panels |
| `02_full_mode_QC.csv` | Per-sample sum of raw weights |
| `03_concordance_full_mode_vs_original.csv` | Agreement between the previous figure values (normalised) and the full-mode rerun |
| `04_group_summary.csv` | Mean, SEM, minimum and maximum proportion per group |
| `04_statistics_DB_vs_ND.csv` | Two-sided Welch's t-test, DB vs ND VCO (n = 3 donors per group) |

## Software

R 4.4.2, Seurat [version], spacexr [version], Matrix, dplyr, tidyr, tibble, ggplot2. Exact versions are written to `results/sessionInfo_*.txt` on each run.

## Licence

MIT (see `LICENSE`).
