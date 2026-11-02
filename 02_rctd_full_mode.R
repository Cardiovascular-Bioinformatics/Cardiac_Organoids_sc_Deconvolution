############################################################################
# 02_rctd_full_mode.R
#
# RCTD (spacexr) deconvolution of VCO bulk RNA-seq in FULL mode.
# Each bulk library is treated as one "pixel". Full mode estimates the
# contribution of every reference cell type.
#
# Inputs : results/Cardiac_organoid_annotated_Seurat.rds   (script 01)
#          data/RKG1vsRKG2.xls  (Novogene read-count table, RK001-RK006)
# Outputs: results/02_full_mode_raw_weights.csv
#          results/02_full_mode_normalised_proportions.csv   <- used in ALL panels
#          results/02_full_mode_QC.csv
#          results/02_rctd_reference.rds, results/02_bulk_spatialRNA.rds
#          results/02_RCTD_full_mode_object.rds
############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
  library(Matrix)
  library(dplyr)
  library(tibble)
})

set.seed(20261102)   # Reference() randomly down-samples to n_max_cells per type
dir.create("results", showWarnings = FALSE)

sample_ids <- sprintf("RK%03d_count", 1:6)   # RK001-003 DB VCO; RK004-006 ND VCO

############################################################################
# 1. Single-cell reference
############################################################################
obj <- readRDS("results/Cardiac_organoid_annotated_Seurat.rds")

sc_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
stopifnot(all(sc_counts@x %% 1 == 0))

cell_types <- obj$celltype
names(cell_types) <- colnames(obj)
keep <- !is.na(cell_types)
sc_counts  <- sc_counts[, keep]
cell_types <- droplevels(factor(cell_types[keep]))

print(table(cell_types))                 # every type must have >= 25 cells
stopifnot(all(table(cell_types) >= 25))

reference <- Reference(counts = sc_counts, cell_types = cell_types,
                       nUMI = colSums(sc_counts), n_max_cells = 10000)

############################################################################
# 2. Bulk RNA-seq counts
############################################################################
bulk <- read.delim("data/RKG1vsRKG2.xls", header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)

bulk_counts <- bulk %>%
  select(gene_name, all_of(sample_ids)) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  group_by(gene_name) %>%                          # collapse duplicated symbols
  summarise(across(everything(), sum), .groups = "drop") %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

stopifnot(all(bulk_counts %% 1 == 0))
print(colSums(bulk_counts))                         # library sizes

bulk_srna <- SpatialRNA(
  coords = data.frame(x = seq_len(ncol(bulk_counts)), y = 1,
                      row.names = colnames(bulk_counts)),
  counts = Matrix(bulk_counts, sparse = TRUE),
  nUMI   = colSums(bulk_counts)
)

saveRDS(reference, "results/02_rctd_reference.rds")
saveRDS(bulk_srna, "results/02_bulk_spatialRNA.rds")

############################################################################
# 3. RCTD, full mode
############################################################################
rctd <- create.RCTD(
  spatialRNA        = bulk_srna,
  reference         = reference,
  max_cores         = 4,
  UMI_min           = 100,
  UMI_max           = 1e9,      # default (2e7) can silently drop deep libraries
  CELL_MIN_INSTANCE = 25        # all other gene-selection parameters: defaults
)
stopifnot(ncol(rctd@spatialRNA@counts) == length(sample_ids))

rctd <- run.RCTD(rctd, doublet_mode = "full")
saveRDS(rctd, "results/02_RCTD_full_mode_object.rds")

############################################################################
# 4. Weights -> proportions
############################################################################
# RCTD weights are not constrained to sum to 1. normalize_weights() divides
# each sample's weights by their sum, giving proportions that sum to 1.
raw_w <- as.matrix(rctd@results$weights)
prop  <- normalize_weights(raw_w)
stopifnot(all(abs(rowSums(prop) - 1) < 1e-8))

conv <- rctd@results$converged
qc <- data.frame(
  sample             = rownames(raw_w),
  sum_of_raw_weights = rowSums(raw_w),
  converged          = if (length(conv) == nrow(raw_w)) as.logical(conv) else NA
)

print(qc)

write.csv(raw_w, "results/02_full_mode_raw_weights.csv")
write.csv(prop,  "results/02_full_mode_normalised_proportions.csv")
write.csv(qc,    "results/02_full_mode_QC.csv", row.names = FALSE)

writeLines(c(capture.output(sessionInfo()),
             paste("spacexr", as.character(packageVersion("spacexr")))),
           "results/sessionInfo_02.txt")
