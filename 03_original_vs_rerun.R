############################################################################
# 03_original_vs_rerun.R
#
# Reproduces the values plotted in the previous version of Fig. 4G-J
# (@results$weights of the original RCTD object, unnormalised) and compares
# them, after normalisation, with the full-mode rerun.
#
# Inputs : data/RKG1_RKG2_RCTD_object.rds  (original RCTD object)
#          results/02_full_mode_normalised_proportions.csv
# Outputs: results/03_original_*.csv
############################################################################

suppressPackageStartupMessages({ library(spacexr); library(Matrix) })

old <- readRDS("data/RKG1_RKG2_RCTD_object.rds")
print(names(old@results))

old_raw <- as.matrix(old@results$weights)
if (is.null(rownames(old_raw))) rownames(old_raw) <- colnames(old@spatialRNA@counts)
old_norm <- normalize_weights(old_raw)

write.csv(old_raw,  "results/03_original_results_weights_raw_used_in_old_Fig4.csv")
write.csv(old_norm, "results/03_original_results_weights_normalised.csv")
write.csv(data.frame(sample = rownames(old_raw), sum_of_raw_weights = rowSums(old_raw)),
          "results/03_original_sum_of_raw_weights.csv", row.names = FALSE)

# Concordance with the full-mode rerun (script 02)
prop <- as.matrix(read.csv("results/02_full_mode_normalised_proportions.csv",
                           row.names = 1, check.names = FALSE))
ct <- intersect(colnames(prop), colnames(old_norm))
concord <- sapply(rownames(prop), function(s)
  cor(prop[s, ct], old_norm[s, ct], method = "pearson"))
print(concord)
write.csv(data.frame(sample = names(concord), pearson_r = concord,
                     max_abs_difference = sapply(rownames(prop), function(s)
                       max(abs(prop[s, ct] - old_norm[s, ct])))),
          "results/03_concordance_full_mode_vs_original.csv", row.names = FALSE)
