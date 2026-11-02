############################################################################
# 01_reference_atlas_GSE263193.R
#
# Re-analysis and annotation of the cardiac organoid scRNA-seq atlas
# (GEO: GSE263193) used as the RCTD reference.
# Produces Supplementary Fig. 5A-C and the annotated reference object.
#
# Inputs  (data/):   GSE263193_barcodes.tsv.gz
#                    GSE263193_genes.tsv.gz
#                    GSE263193_seurat_sf_counts.mtx.gz
#                    GSE263193_seurat_metadata.tsv.gz
# Outputs (results/): Cardiac_organoid_annotated_Seurat.rds
#                    01_reference_cell_annotations.csv  (barcode -> cell type)
#                    01_reference_cells_per_celltype.csv
#                    01_all_cluster_markers.csv, 01_top10_cluster_markers.csv
#                    SuppFig5A_UMAP.png, SuppFig5B_known_markers.png,
#                    SuppFig5C_top10_markers.png
#
# Run from the repository root:  Rscript R/01_reference_atlas_GSE263193.R
############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})

set.seed(20261102)
dir.create("results", showWarnings = FALSE)

# Optional: cluster-level SingleR (HPCA) labels and lineage module scores
# were used as supporting evidence for the manual annotation below.
run_supporting_annotation <- FALSE

############################################################################
# 1. Load GSE263193
############################################################################
barcodes <- read.delim(gzfile("data/GSE263193_barcodes.tsv.gz"),
                       header = FALSE, stringsAsFactors = FALSE)[, 1]
genes    <- read.delim(gzfile("data/GSE263193_genes.tsv.gz"),
                       header = FALSE, stringsAsFactors = FALSE)
gene_names <- if (ncol(genes) == 1) genes[, 1] else genes[, 2]

counts <- readMM(gzfile("data/GSE263193_seurat_sf_counts.mtx.gz"))
rownames(counts) <- make.unique(gene_names)
colnames(counts) <- barcodes

# RCTD requires raw integer counts
stopifnot(all(counts@x %% 1 == 0))

obj <- CreateSeuratObject(counts = counts, project = "Cardiac_Organoid",
                          min.cells = 3, min.features = 200)

metadata <- read.delim(gzfile("data/GSE263193_seurat_metadata.tsv.gz"),
                       row.names = 1, stringsAsFactors = FALSE)
obj <- AddMetaData(obj, metadata = metadata[colnames(obj), ])

############################################################################
# 2. Quality control
############################################################################
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 &
                            percent.mt < 15)

############################################################################
# 3. Normalisation, dimensionality reduction, clustering
############################################################################
obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, features = VariableFeatures(obj))
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.2)          # Louvain (Seurat default)
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca")

############################################################################
# 4. Optional supporting annotation (SingleR / module scores)
############################################################################
if (run_supporting_annotation) {
  suppressPackageStartupMessages({ library(SingleR); library(celldex) })
  hpca <- HumanPrimaryCellAtlasData()
  pred <- SingleR(test = LayerData(obj, assay = "RNA", layer = "data"),
                  ref = hpca, labels = hpca$label.main,
                  clusters = obj$RNA_snn_res.0.2)
  write.csv(data.frame(cluster = rownames(pred), HPCA_label = pred$labels,
                       HPCA_pruned = pred$pruned.labels,
                       delta = pred$delta.next),
            "results/01_SingleR_HPCA_cluster_labels.csv", row.names = FALSE)
}

############################################################################
# 5. Manual annotation of resolution-0.2 clusters
############################################################################
cluster_annotations <- c(
  "1"  = "Cardiomyocytes_1",              # NPPA/NPPB stress-responsive signature
  "6"  = "Cardiomyocytes_2",
  "0"  = "Developing_Cardiomyocytes",
  "8"  = "Mural_Cells",
  "2"  = "Fibroblasts1",
  "3"  = "Fibroblasts2",
  "11" = "Fibroblast_Progenitors",
  "4"  = "Embryonic_Cardiac_Progenitors1",
  "5"  = "Embryonic_Cardiac_Progenitors2",
  "7"  = "Embryonic_Stem_Cells_1",
  "9"  = "Embryonic_Stem_Cells_2",
  "10" = "Endothelial_cells"
)

celltype_order <- c(
  "Cardiomyocytes_1", "Cardiomyocytes_2", "Developing_Cardiomyocytes",
  "Mural_Cells", "Fibroblasts1", "Fibroblasts2", "Fibroblast_Progenitors",
  "Embryonic_Cardiac_Progenitors1", "Embryonic_Cardiac_Progenitors2",
  "Embryonic_Stem_Cells_1", "Embryonic_Stem_Cells_2", "Endothelial_cells"
)

obj$celltype <- factor(unname(cluster_annotations[as.character(obj$RNA_snn_res.0.2)]),
                       levels = celltype_order)
stopifnot(!any(is.na(obj$celltype)))
Idents(obj) <- obj$celltype

############################################################################
# 6. Supplementary Fig. 5A-C
############################################################################
p_umap <- DimPlot(obj, reduction = "umap", group.by = "celltype",
                  label = TRUE, repel = TRUE, label.size = 4) +
  ggtitle("Cardiac organoid cell atlas (GSE263193)") + theme_classic()
ggsave("results/SuppFig5A_UMAP.png", p_umap, width = 10, height = 8, dpi = 600)

known_markers <- c(
  "TNNT2", "ACTN2", "MYH6", "MYH7", "TNNI3", "RYR2", "PLN",     # cardiomyocytes
  "MYL4", "TNNI1", "NPPA", "NPPB", "ISL1",                      # developing / stress
  "PDGFRB", "RGS5", "CSPG4", "MCAM", "DES", "ACTA2",            # mural
  "PDGFRA", "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "TCF21",# fibroblast
  "NKX2-5", "MEF2C", "HAND1", "TBX5",                           # cardiac progenitor
  "POU5F1", "NANOG", "SOX2", "LIN28A", "KLF4",                  # stem-like
  "PECAM1", "VWF", "KDR", "EMCN", "ENG"                         # endothelial
)
known_markers <- unique(known_markers[known_markers %in% rownames(obj)])

p_known <- DotPlot(obj, features = known_markers, group.by = "celltype",
                   dot.scale = 8) +
  RotatedAxis() + ggtitle("Known lineage marker expression") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))
ggsave("results/SuppFig5B_known_markers.png", p_known,
       width = 20, height = 10, dpi = 600)

markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10)

write.csv(markers, "results/01_all_cluster_markers.csv", row.names = FALSE)
write.csv(top10,   "results/01_top10_cluster_markers.csv", row.names = FALSE)

p_top10 <- DotPlot(obj, features = unique(top10$gene), group.by = "celltype",
                   dot.scale = 6) +
  ggtitle("Top 10 markers per cluster") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
ggsave("results/SuppFig5C_top10_markers.png", p_top10,
       width = 24, height = 12, dpi = 600, limitsize = FALSE)

############################################################################
# 7. Save annotated reference
############################################################################
# The barcode -> cell type table allows the reference to be rebuilt exactly,
# independent of cluster numbering in other Seurat versions.
write.csv(data.frame(barcode = colnames(obj),
                     cluster_res0.2 = obj$RNA_snn_res.0.2,
                     celltype = obj$celltype),
          "results/01_reference_cell_annotations.csv", row.names = FALSE)
write.csv(as.data.frame(table(celltype = obj$celltype)),
          "results/01_reference_cells_per_celltype.csv", row.names = FALSE)

saveRDS(obj, "results/Cardiac_organoid_annotated_Seurat.rds")
writeLines(capture.output(sessionInfo()), "results/sessionInfo_01.txt")
