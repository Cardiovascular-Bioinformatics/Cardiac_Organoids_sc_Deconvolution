############################################################################
# 04_figures.R
#
# All Fig. 4G-J panels are drawn from ONE matrix: the normalised full-mode RCTD
# proportions (rows sum to 1) from script 02.
#
# Outputs: results/Fig4G_composition.png/.pdf
#          results/Fig4H_Endothelial_cells.png/.pdf
#          results/Fig4I_Cardiomyocytes_1.png/.pdf
#          results/Fig4J_Cardiomyocytes_2.png/.pdf
#          results/04_group_summary.csv
#          results/04_statistics_DB_vs_ND.csv
############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

sample_info <- data.frame(
  sample = sprintf("RK%03d_count", 1:6),
  donor  = sprintf("RK%03d", 1:6),
  group  = factor(rep(c("DB VCO", "ND VCO"), each = 3),
                  levels = c("ND VCO", "DB VCO"))
)
group_cols <- c("ND VCO" = "#F7A8B8", "DB VCO" = "#F7941D")  # as in Fig. 4A

# Same order as obj$celltype in script 01, so colours match Supp. Fig. 5A
celltype_order <- c(
  "Cardiomyocytes_1", "Cardiomyocytes_2", "Developing_Cardiomyocytes",
  "Mural_Cells", "Fibroblasts1", "Fibroblasts2", "Fibroblast_Progenitors",
  "Embryonic_Cardiac_Progenitors1", "Embryonic_Cardiac_Progenitors2",
  "Embryonic_Stem_Cells_1", "Embryonic_Stem_Cells_2", "Endothelial_cells"
)

prop  <- as.matrix(read.csv("results/02_full_mode_normalised_proportions.csv",
                            row.names = 1, check.names = FALSE))
raw_w <- as.matrix(read.csv("results/02_full_mode_raw_weights.csv",
                            row.names = 1, check.names = FALSE))

save_both <- function(name, plot, width, height) {
  ggsave(file.path("results", paste0(name, ".png")), plot,
         width = width, height = height, dpi = 600)
  ggsave(file.path("results", paste0(name, ".pdf")), plot,
         width = width, height = height)
}

############################################################################
# 1. Long format, total cardiomyocytes, group summaries and Welch t-tests
############################################################################
prop_long <- as.data.frame(prop) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "celltype", values_to = "proportion") %>%
  left_join(sample_info, by = "sample") %>%
  mutate(celltype = factor(celltype, levels = celltype_order))

composites <- prop_long %>%
  group_by(sample, donor, group) %>%
  summarise(
    All_Cardiomyocytes = sum(proportion[celltype %in%
                           c("Cardiomyocytes_1", "Cardiomyocytes_2")]),
    .groups = "drop") %>%
  pivot_longer(All_Cardiomyocytes, names_to = "celltype",
               values_to = "proportion")

all_long <- bind_rows(mutate(prop_long, celltype = as.character(celltype)),
                      composites)

group_summary <- all_long %>%
  group_by(celltype, group) %>%
  summarise(mean = mean(proportion), sem = sd(proportion) / sqrt(n()),
            min = min(proportion), max = max(proportion), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = c(mean, sem, min, max),
              names_glue = "{group}_{.value}")
write.csv(group_summary, "results/04_group_summary.csv", row.names = FALSE)

# Two-sided Welch's t-test, DB vs ND VCO (independent donors, n = 3 per group)
welch_p <- function(y, g) tryCatch(t.test(y ~ g, var.equal = FALSE)$p.value,
                                   error = function(e) NA_real_)
statistics <- all_long %>%
  group_by(celltype) %>%
  summarise(mean_ND = mean(proportion[group == "ND VCO"]),
            mean_DB = mean(proportion[group == "DB VCO"]),
            p_welch = welch_p(proportion, group), .groups = "drop") %>%
  arrange(p_welch)
print(statistics, n = Inf)
write.csv(statistics, "results/04_statistics_DB_vs_ND.csv", row.names = FALSE)

############################################################################
# 2. Fig. 4G: per-donor composition
############################################################################
fig_g <- ggplot(prop_long, aes(donor, proportion, fill = celltype)) +
  geom_col(position = "stack", width = 0.85) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "RCTD-estimated proportion", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12))
save_both("Fig4G_composition", fig_g, 8, 5)

############################################################################
# 3. Fig. 4H-J: per-population boxplots
############################################################################
make_box <- function(ct) {
  df <- filter(all_long, celltype == ct)
  ggplot(df, aes(group, proportion, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.5) +
    geom_point(size = 2.5, position = position_jitter(width = 0.08, seed = 1)) +
    scale_fill_manual(values = group_cols) +
    labs(title = ct, x = NULL, y = "RCTD-estimated proportion") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          axis.text = element_text(colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
}

main_panels <- c(H = "Endothelial_cells", I = "Cardiomyocytes_1",
                 J = "Cardiomyocytes_2")
for (k in names(main_panels)) {
  save_both(sprintf("Fig4%s_%s", k, main_panels[[k]]),
            make_box(main_panels[[k]]), 3.2, 4.2)
}

writeLines(capture.output(sessionInfo()), "results/sessionInfo_04.txt")
