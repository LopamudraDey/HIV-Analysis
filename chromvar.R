library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)
library(chromVAR)
library(motifmatchr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# ---------------------------------------------------------------
# Settings
# ---------------------------------------------------------------
obj_path      <- "E:/multiome_with_Marie_celltype_atac.rds"
out_dir       <- "E:/HIV vs Mock_ Clusterwise_MotifEnrichment/Heatmap1"
cluster_col   <- "seurat_clusters"   # or your cell-type column, e.g. "Marie_celltype"
condition_col <- "condition"
hiv_label     <- "HIV"
mock_label    <- "Mock"

# Columns of the heatmap:
#   "cluster"           -> one column per cluster (HIV + ssDNA pooled)
#   "cluster_condition" -> one column per cluster x condition
column_mode   <- "cluster"

n_top_per_col <- 8      # top motifs picked per column
min_cells     <- 10     # min cells for a column to be kept
color_lim     <- 1.5    # colour scale limit (row z-score)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------
# 1. Load object, add motifs, run chromVAR (once)
# ---------------------------------------------------------------
obj <- readRDS(obj_path)
DefaultAssay(obj) <- "ATAC"

if (is.null(Motifs(obj))) {
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
  )
  obj <- AddMotifs(obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
}

if (!"chromvar" %in% names(obj@assays)) {
  message("Running chromVAR manually (this can take a while)...")
  
  library(SummarizedExperiment)
  library(Matrix)
  library(BiocParallel)
  register(SerialParam())   # multicore is not supported on Windows
  
  # 1. Peak counts + peak ranges -> SummarizedExperiment
  counts_mat  <- GetAssayData(obj, assay = "ATAC", layer = "counts")
  peak_ranges <- granges(obj[["ATAC"]])
  names(peak_ranges) <- rownames(counts_mat)
  
  se <- SummarizedExperiment(
    assays    = list(counts = counts_mat),
    rowRanges = peak_ranges,
    colData   = obj[[]]
  )
  
  # drop peaks with no counts (chromVAR errors / warns on these)
  se <- se[Matrix::rowSums(assay(se, "counts")) > 0, ]
  
  # 2. Motif matches
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
  )
  motif_ix <- matchMotifs(pfm, se, genome = BSgenome.Hsapiens.UCSC.hg38)
  
  # 3. chromVAR
  se  <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)
  bg  <- getBackgroundPeaks(se)
  dev <- computeDeviations(object = se, annotations = motif_ix, background_peaks = bg)
  
  # 4. z-scores (note: deviationScores, not deviationsScores)
  dev_scores <- deviationScores(dev)
  dev_scores[is.na(dev_scores)] <- 0
  
  # rownames are JASPAR IDs -> convert to TF names
  rownames(dev_scores) <- make.unique(TFBSTools::name(pfm)[rownames(dev_scores)])
  
  # 5. Add back to the Seurat object
  obj[["chromvar"]] <- CreateAssayObject(data = dev_scores)
  saveRDS(obj, file.path(out_dir, "obj_with_chromvar.rds"))
  manual_names <- TRUE
}

# ---------------------------------------------------------------
# 2. Average chromVAR deviation z-scores per column group
# ---------------------------------------------------------------
md <- obj[[]]
keep <- as.character(md[[condition_col]]) %in% c(hiv_label, mock_label)

grp <- if (column_mode == "cluster") {
  as.character(md[[cluster_col]])
} else {
  paste0(as.character(md[[cluster_col]]), "_", as.character(md[[condition_col]]))
}
names(grp) <- rownames(md)
grp <- grp[keep]

mat <- GetAssayData(obj, assay = "chromvar", layer = "data")[, names(grp)]

groups <- names(which(table(grp) >= min_cells))
# natural ordering (numeric clusters first)
groups <- groups[order(suppressWarnings(as.numeric(sub("_.*", "", groups))), groups)]

avg <- sapply(groups, function(g) Matrix::rowMeans(mat[, grp == g, drop = FALSE]))

# Motif IDs -> TF names
if (all(grepl("^MA[0-9]+", rownames(avg)))) {   # still JASPAR IDs (e.g. from RunChromVAR)
  rownames(avg) <- make.unique(ConvertMotifID(obj, id = rownames(avg), assay = "ATAC"))
}

write.csv(avg, file.path(out_dir, "chromVAR_mean_deviation_per_group.csv"))

# ---------------------------------------------------------------
# 3. Row z-score across columns + pick top motifs per column
# ---------------------------------------------------------------
z <- t(scale(t(avg)))
z <- z[complete.cases(z), , drop = FALSE]

picked <- unlist(lapply(colnames(z), function(g) {
  names(sort(z[, g], decreasing = TRUE))[seq_len(n_top_per_col)]
}))
picked <- unique(picked)

# order rows by the column where they peak, then by peak value
peak_col <- apply(z[picked, , drop = FALSE], 1, which.max)
peak_val <- apply(z[picked, , drop = FALSE], 1, max)
picked   <- picked[order(peak_col, -peak_val)]
zp       <- z[picked, , drop = FALSE]

# ---------------------------------------------------------------
# 4. Draw heatmap (blue - yellow - red, like the reference panel)
# ---------------------------------------------------------------
col_fun <- colorRamp2(c(-color_lim, 0, color_lim),
                      c("#4575B4", "#FFFFBF", "#D73027"))

ht <- Heatmap(
  zp,
  name = "Row z-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 90,
  column_title = if (column_mode == "cluster") "Motif activity (chromVAR) by cluster"
  else "Motif activity (chromVAR) by cluster x condition",
  column_title_gp = gpar(fontface = "bold", fontsize = 12),
  border = FALSE,
  rect_gp = gpar(col = NA),
  heatmap_legend_param = list(at = c(-color_lim, 0, color_lim))
)

h <- max(5, 0.16 * nrow(zp) + 2)
w <- max(5, 0.45 * ncol(zp) + 3)

pdf(file.path(out_dir, "Motif_activity_heatmap.pdf"), width = w, height = h)
draw(ht)
dev.off()

png(file.path(out_dir, "Motif_activity_heatmap.png"), width = w, height = h,
    units = "in", res = 300)
draw(ht)
dev.off()

message("Done. Heatmap saved in: ", out_dir)




library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)

# ---------------------------------------------------------------
# Settings
# ---------------------------------------------------------------
out_dir       <- "E:/HIV vs Mock_ Clusterwise_MotifEnrichment/Heatmap"
obj_path      <- file.path(out_dir, "obj_with_chromvar.rds")   # saved by the previous script
cluster_col   <- "seurat_clusters"
condition_col <- "condition"
hiv_label     <- "HIV"
mock_label    <- "ssDNA"

n_top         <- 6      # motifs shown per panel in EACH direction (higher in HIV / higher in ssDNA)
min_cells     <- 10     # min cells per condition within a cluster
padj_cut      <- 0.05   # motifs marked with * if adjusted P (Wilcoxon, HIV vs ssDNA) below this
panel_ncol    <- 4      # number of panels per row in the grid

# ---------------------------------------------------------------
# 1. Load stored object + chromVAR matrix
# ---------------------------------------------------------------
obj <- readRDS(obj_path)
mat <- GetAssayData(obj, assay = "chromvar", layer = "data")

# make sure rows are TF names (not JASPAR IDs)
if (all(grepl("^MA[0-9]+", rownames(mat)))) {
  rownames(mat) <- make.unique(ConvertMotifID(obj, id = rownames(mat), assay = "ATAC"))
}

md <- obj[[]]
md <- md[as.character(md[[condition_col]]) %in% c(hiv_label, mock_label), , drop = FALSE]
cl_ids <- sort(unique(as.character(md[[cluster_col]])))
cl_ids <- cl_ids[order(suppressWarnings(as.numeric(cl_ids)), cl_ids)]

# ---------------------------------------------------------------
# 2. Per cluster: mean chromVAR z-score per condition, delta, Wilcoxon test
# ---------------------------------------------------------------
stats_list <- list()
for (cl in cl_ids) {
  cells_h <- rownames(md)[as.character(md[[cluster_col]]) == cl & as.character(md[[condition_col]]) == hiv_label]
  cells_m <- rownames(md)[as.character(md[[cluster_col]]) == cl & as.character(md[[condition_col]]) == mock_label]
  if (length(cells_h) < min_cells || length(cells_m) < min_cells) {
    message("Skipping cluster ", cl, " (HIV: ", length(cells_h), ", ssDNA: ", length(cells_m), " cells)")
    next
  }
  x <- as.matrix(mat[, cells_h, drop = FALSE])
  y <- as.matrix(mat[, cells_m, drop = FALSE])
  
  p <- vapply(seq_len(nrow(x)), function(i) {
    tryCatch(wilcox.test(x[i, ], y[i, ])$p.value, error = function(e) NA_real_)
  }, numeric(1))
  
  stats_list[[cl]] <- data.frame(
    cluster = cl,
    motif   = rownames(x),
    mean_HIV   = rowMeans(x),
    mean_ssDNA = rowMeans(y),
    delta      = rowMeans(x) - rowMeans(y),
    p_val      = p,
    p_adj      = p.adjust(p, method = "BH"),
    stringsAsFactors = FALSE
  )
}
stats <- bind_rows(stats_list)
write.csv(stats, file.path(out_dir, "chromVAR_HIV_vs_ssDNA_stats_per_cluster.csv"), row.names = FALSE)

# ---------------------------------------------------------------
# 3. Pick top motifs per cluster panel
# ---------------------------------------------------------------
top <- stats %>%
  filter(!is.na(delta)) %>%
  group_by(cluster) %>%
  arrange(desc(delta), .by_group = TRUE) %>%
  filter(row_number() <= n_top | row_number() > n() - n_top) %>%
  mutate(rank = row_number()) %>%      # 1 = most HIV-high ... last = most ssDNA-high
  ungroup()

plot_df <- top %>%
  tidyr::pivot_longer(c(mean_HIV, mean_ssDNA), names_to = "condition", values_to = "value") %>%
  mutate(
    condition = factor(sub("mean_", "", condition), levels = c(mock_label, hiv_label)),
    key       = paste(cluster, motif, sep = "__"),
    label     = paste0(motif, ifelse(!is.na(p_adj) & p_adj < padj_cut, " *", "")),
    cluster_f = factor(paste("Cluster", cluster), levels = paste("Cluster", cl_ids))
  )

# row order inside each panel: HIV-high on top
key_levels <- plot_df %>% distinct(cluster, rank, key) %>% arrange(cluster, desc(rank)) %>% pull(key)
plot_df$key <- factor(plot_df$key, levels = unique(key_levels))
lab_map <- setNames(plot_df$label, as.character(plot_df$key))

lim <- round(as.numeric(quantile(abs(plot_df$value), 0.98)), 1)

# ---------------------------------------------------------------
# 4. Multi-panel heatmap grid
# ---------------------------------------------------------------
p_grid <- ggplot(plot_df, aes(condition, key, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  facet_wrap(~cluster_f, scales = "free_y", ncol = panel_ncol) +
  scale_y_discrete(labels = lab_map, position = "right") +
  scale_fill_gradientn(
    colours = c("#4575B4", "#FFFFBF", "#D73027"),
    limits = c(-lim, lim), oob = scales::squish,
    name = "Mean chromVAR\nz-score"
  ) +
  labs(
    x = NULL, y = NULL,
    title = paste0("Top differential motif activity per cluster: ", hiv_label, " vs ", mock_label),
    subtitle = paste0("Top ", n_top, " higher in ", hiv_label, " (upper) and ", n_top, " higher in ",
                      mock_label, " (lower) per cluster;  * = adj. P < ", padj_cut, " (Wilcoxon, BH)")
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid  = element_blank(),
    axis.text.y = element_text(size = 6.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    strip.text  = element_text(face = "bold", size = 11),
    plot.title  = element_text(face = "bold"),
    legend.position = "right"
  )

n_panels <- length(unique(plot_df$cluster))
n_rows   <- ceiling(n_panels / panel_ncol)
w <- 3.0 * min(panel_ncol, n_panels) + 1.5
h <- n_rows * (0.17 * 2 * n_top + 1.3) + 1

ggsave(file.path(out_dir, "chromVAR_cluster_subplots.pdf"), p_grid, width = w, height = h, limitsize = FALSE)
ggsave(file.path(out_dir, "chromVAR_cluster_subplots.png"), p_grid, width = w, height = h, dpi = 300, limitsize = FALSE)

message("Done. Saved subplots and stats table in: ", out_dir)



##################################################################################################
# chromvar Cluster wise ###############################################################################################
out_dir       <- "E:/MotifEnrichment2/HIV_vs_ssDNA_chromVAR"
cluster_col   <- "seurat_clusters"
cluster_id    <- "0"
condition_col <- "condition"
hiv_label     <- "HIV"
mock_label    <- "ssDNA"

min_cells_per_peak <- 10     # keep peaks detected in at least this many cluster-0 cells
padj_cut      <- 0.05
min_delta     <- 0.3         # minimum |HIV - ssDNA| difference in mean chromVAR z-score
n_label       <- 12          # labelled motifs per direction

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
register(SerialParam())      # multicore is not supported on Windows

# ---------------------------------------------------------------
# 1. Subset cluster 0 (HIV + ssDNA cells only)
# ---------------------------------------------------------------
obj <- readRDS(obj_path)
DefaultAssay(obj) <- "ATAC"

md <- obj[[]]
cells_keep <- rownames(md)[
  as.character(md[[cluster_col]]) == cluster_id &
    as.character(md[[condition_col]]) %in% c(hiv_label, mock_label)
]
sub_obj <- subset(obj, cells = cells_keep)
rm(obj); gc()

message("Cluster ", cluster_id, " cells per condition:")
print(table(sub_obj[[condition_col]][, 1]))

# ---------------------------------------------------------------
# 2. Build SummarizedExperiment from cluster-0 peak counts only
# ---------------------------------------------------------------
counts_mat  <- GetAssayData(sub_obj, assay = "ATAC", layer = "counts")
peak_ranges <- granges(sub_obj[["ATAC"]])
names(peak_ranges) <- rownames(counts_mat)

# drop empty cells and rarely detected peaks (re-computed for this cluster only)
counts_mat <- counts_mat[, Matrix::colSums(counts_mat) > 0, drop = FALSE]
keep_peaks <- Matrix::rowSums(counts_mat > 0) >= min_cells_per_peak
counts_mat <- counts_mat[keep_peaks, ]
peak_ranges <- peak_ranges[rownames(counts_mat)]
message("Using ", nrow(counts_mat), " peaks and ", ncol(counts_mat), " cells")

se <- SummarizedExperiment(
  assays    = list(counts = counts_mat),
  rowRanges = peak_ranges,
  colData   = sub_obj[[]][colnames(counts_mat), , drop = FALSE]
)

# ---------------------------------------------------------------
# 3. chromVAR from scratch on this cluster
# ---------------------------------------------------------------
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
)
motif_ix <- matchMotifs(pfm, se, genome = BSgenome.Hsapiens.UCSC.hg38)

se  <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)
bg  <- getBackgroundPeaks(se)
dev <- computeDeviations(object = se, annotations = motif_ix, background_peaks = bg)

dev_scores <- deviationScores(dev)
dev_scores[is.na(dev_scores)] <- 0

# JASPAR IDs -> TF names
id2name <- setNames(TFBSTools::name(pfm), names(pfm))
ids <- rownames(dev_scores)
rownames(dev_scores) <- make.unique(ifelse(ids %in% names(id2name), id2name[ids], ids))

#saveRDS(dev_scores, file.path(out_dir, "cluster1_chromVAR_zscores.rds"))

# ---------------------------------------------------------------
# 4. HIV vs ssDNA within cluster 0
# ---------------------------------------------------------------
cond <- as.character(se$condition)
if (is.null(se$condition)) cond <- as.character(colData(se)[[condition_col]])
cells_h <- colnames(dev_scores)[cond == hiv_label]
cells_m <- colnames(dev_scores)[cond == mock_label]

x <- dev_scores[, cells_h, drop = FALSE]
y <- dev_scores[, cells_m, drop = FALSE]

p <- vapply(seq_len(nrow(x)), function(i) {
  tryCatch(wilcox.test(x[i, ], y[i, ])$p.value, error = function(e) NA_real_)
}, numeric(1))

res <- data.frame(
  motif      = rownames(x),
  mean_HIV   = rowMeans(x),
  mean_ssDNA = rowMeans(y),
  delta      = rowMeans(x) - rowMeans(y),
  p_val      = p,
  p_adj      = p.adjust(p, method = "BH"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    significance = case_when(
      p_adj < padj_cut & delta >=  min_delta ~ paste("Higher in", hiv_label),
      p_adj < padj_cut & delta <= -min_delta ~ paste("Higher in", mock_label),
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(desc(delta))

write.csv(res, file.path(out_dir, "cluster0_chromVAR_HIV_vs_ssDNA.csv"), row.names = FALSE)
message("Cells: ", length(cells_h), " HIV, ", length(cells_m), " ssDNA")
print(table(res$significance))

# ---------------------------------------------------------------
# 5. Plots
# ---------------------------------------------------------------
cols <- c("#C43C39", "#3676A8", "grey75")
names(cols) <- c(paste("Higher in", hiv_label), paste("Higher in", mock_label), "Not significant")

# 5a. volcano: effect size vs significance
lab <- res %>%
  filter(significance != "Not significant") %>%
  group_by(significance) %>%
  slice_max(order_by = abs(delta), n = n_label, with_ties = FALSE) %>%
  ungroup()

p_vol <- ggplot(res, aes(delta, -log10(pmax(p_adj, 1e-300)), colour = significance)) +
  geom_point(size = 1.2, alpha = 0.8) +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = c(-min_delta, min_delta), linetype = "dotted", linewidth = 0.4) +
  geom_text_repel(data = lab, aes(label = motif), size = 3, max.overlaps = Inf,
                  box.padding = 0.35, show.legend = FALSE) +
  scale_colour_manual(values = cols) +
  labs(
    title = paste0("Cluster ", cluster_id, ": motif activity, ", hiv_label, " vs ", mock_label),
    subtitle = "chromVAR computed on cluster 0 cells only",
    x = paste0("Difference in mean chromVAR z-score (", hiv_label, " - ", mock_label, ")"),
    y = expression(-log[10]("adj. P")), colour = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave(file.path(out_dir, "cluster0_volcano.png"), p_vol, width = 6.5, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "cluster0_volcano.pdf"), p_vol, width = 6.5, height = 5.5)

# 5b. bar plot of the top motifs in each direction
bars <- res %>%
  filter(significance != "Not significant") %>%
  group_by(significance) %>%
  slice_max(order_by = abs(delta), n = n_label, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(motif = factor(motif, levels = motif[order(delta)]))

if (nrow(bars) > 0) {
  p_bar <- ggplot(bars, aes(delta, motif, fill = significance)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = cols, name = NULL) +
    labs(
      title = paste0("Cluster ", cluster_id, ": top differential motifs"),
      x = paste0("Difference in mean chromVAR z-score (", hiv_label, " - ", mock_label, ")"),
      y = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
  
  h <- max(4, 0.25 * nrow(bars) + 1.8)
  ggsave(file.path(out_dir, "cluster0_top_motifs_bar.png"), p_bar, width = 6, height = h, dpi = 300)
  ggsave(file.path(out_dir, "cluster0_top_motifs_bar.pdf"), p_bar, width = 6, height = h)
}

message("Done. Results saved in: ", out_dir)

library(tidyr)

# Select top 6 significant motifs (3 HIV-enriched, 3 Mock-enriched)
top_6_motifs <- res %>%
  filter(significance != "Not significant") %>%
  group_by(significance) %>%
  slice_min(order_by = p_adj, n = 10) %>%
  pull(motif)

# Reshape dev_scores for plotting
df_violin <- as.data.frame(t(dev_scores[top_6_motifs, , drop = FALSE]))
df_violin$cell <- rownames(df_violin)
df_violin$condition <- cond
df_long <- pivot_longer(df_violin, cols = all_of(top_6_motifs), names_to = "motif", values_to = "score")

p_violin <- ggplot(df_long, aes(x = condition, y = score, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, alpha = 0.8) +
  facet_wrap(~motif, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("#C43C39", "#3676A8")) +
  labs(
    title = paste0("Cluster ", cluster_id, ": Top Motif Activity Distribution"),
    x = NULL, y = "chromVAR Z-score", fill = "Condition"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

p_violin

ggsave(file.path(out_dir, "cluster0_top_motifs_violin.png"), p_violin, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "cluster0_top_motifs_violin.pdf"), p_violin, width = 8, height = 6)



############################Cluster wise Plots################################################################



library(ggplot2)
library(dplyr)
library(tidyr)
library(Seurat)

# ---------------------------------------------------------------
# 1. Load Data & Extract Metadata
# ---------------------------------------------------------------
obj <- readRDS("E:/multiome_with_Marie_celltype_atac.rds")

cluster_col   <- "seurat_clusters"
condition_col <- "condition" # Ensure this matches your metadata column name

# Extract metadata table
md <- obj[[]]

# ---------------------------------------------------------------
# 2. Summarize Cell Counts per Cluster & Condition
# ---------------------------------------------------------------
cell_counts <- md %>%
  group_by(.data[[cluster_col]], .data[[condition_col]]) %>%
  summarise(count = n(), .groups = "drop") %>%
  rename(cluster = all_of(cluster_col), condition = all_of(condition_col))

# Ensure all clusters and conditions are factorized properly for plotting
cell_counts$cluster <- factor(cell_counts$cluster, levels = sort(unique(as.numeric(as.character(cell_counts$cluster)))))

# Define custom color palette for HIV, Mock, and ssDNA
condition_colors <- c(
  "HIV"   = "#C43C39",
  "Mock"  = "#3676A8",
  "ssDNA" = "#E69F00"
)

# ---------------------------------------------------------------
# 3. Option A: Dodged Bar Plot (All clusters side-by-side)
# ---------------------------------------------------------------
p_dodge <- ggplot(cell_counts, aes(x = cluster, y = count, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single"), width = 0.7) +
  scale_fill_manual(values = condition_colors) +
  labs(
    title = "Cell Counts per Cluster Across Conditions",
    x = "Seurat Cluster",
    y = "Number of Cells",
    fill = "Condition"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "right"
  )
p_dodge
p_boxed_facets <- ggplot(cell_counts, aes(x = condition, y = count, fill = condition)) +
  geom_bar(stat = "identity", width = 0.65, color = "black", linewidth = 0.3) +
  facet_wrap(~ cluster, scales = "free_y", nrow = 3) +
  scale_fill_manual(values = condition_colors) +
  labs(
    title = "Cell Counts per Condition Across All 19 Clusters",
    x = NULL,
    y = "Number of Cells",
    fill = "Condition"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    
    # --- Boxed Cluster Headings (Strip Labels) ---
    strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 10, color = "black"),
    
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )
ggsave("E:/MotifEnrichment/Cluster_Condition_Counts_Dodged.png", p_dodge, width = 10, height = 5, dpi = 300)
ggsave("E:/MotifEnrichment/Cluster_Condition_Counts_Dodged.pdf", p_dodge, width = 10, height = 5)
p_boxed_facets
# ---------------------------------------------------------------
# 4. Option B: Faceted Panel Plot (Similar to your reference style)
# ---------------------------------------------------------------
p_facet <- ggplot(cell_counts, aes(x = condition, y = count, fill = condition)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_wrap(~ cluster, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = condition_colors) +
  labs(
    title = "Cell Counts per Condition by Cluster",
    x = NULL,
    y = "Number of Cells",
    fill = "Condition"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )
p_facet
if (!"nCount_ATAC" %in% colnames(md)) {
  stop("nCount_ATAC not found in Seurat metadata. Ensure your ATAC assay is processed.")
}

# ---------------------------------------------------------------
# 2. Summarize ATAC-seq Counts per Cluster & Condition
# ---------------------------------------------------------------
atac_summary <- md %>%
  group_by(.data[[cluster_col]], .data[[condition_col]]) %>%
  summarise(
    total_atac_counts = sum(nCount_ATAC, na.rm = TRUE), # Total fragments
    mean_atac_counts  = mean(nCount_ATAC, na.rm = TRUE), # Average depth per cell
    cell_count        = n(),
    .groups = "drop"
  ) %>%
  rename(cluster = all_of(cluster_col), condition = all_of(condition_col))

# Order clusters numerically
atac_summary$cluster <- factor(atac_summary$cluster, levels = sort(unique(as.numeric(as.character(atac_summary$cluster)))))

# Define custom color palette for HIV, Mock, and ssDNA
condition_colors <- c(
  "HIV"   = "#C43C39",
  "Mock"  = "#3676A8",
  "ssDNA" = "#E69F00"
)

# ---------------------------------------------------------------
# 3. Generate Facet Plot for Total ATAC Fragment Counts
# ---------------------------------------------------------------
p_atac_total <- ggplot(atac_summary, aes(x = condition, y = total_atac_counts / 1e3, fill = condition)) +
  geom_bar(stat = "identity", width = 0.65, color = "black", linewidth = 0.3) +
  facet_wrap(~ cluster, scales = "free_y", nrow = 3) +
  scale_fill_manual(values = condition_colors) +
  labs(
    title = "Total ATAC-seq Fragment Counts (Thousands) per Condition Across Clusters",
    x = NULL,
    y = "Total ATAC Fragments (in thousands)",
    fill = "Condition"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    
    # Boxed Cluster Headings
    strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 10, color = "black"),
    
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )
p_atac_total
ggsave("E:/MotifEnrichment/Cluster_Condition_Counts_Faceted.png", p_facet, width = 12, height = 6, dpi = 300)
ggsave("E:/MotifEnrichment/Cluster_Condition_Counts_Faceted.pdf", p_facet, width = 12, height = 6)

message("Cell count plots generated and saved successfully!")



library(Seurat)
library(Signac)
library(chromVAR)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(GenomicRanges)
library(BiocParallel)
library(dplyr)
library(ggplot2)
library(ggrepel)

# ---------------- Settings ----------------
obj_path      <- "E:/multiome_with_Marie_celltype_atac.rds"
out_dir       <- "E:/MotifEnrichment2/HIV_vs_Mock_chromVAR"
cluster_col   <- "seurat_clusters"
condition_col <- "condition"
hiv_label     <- "HIV"
mock_label    <- "Mock"

min_cells_per_peak      <- 10   # keep peaks detected in at least this many cells of the cluster
min_cells_per_condition <- 10   # skip cluster if either condition has fewer cells
padj_cut  <- 0.05
min_delta <- 0.3                # minimum |HIV - ssDNA| difference in mean chromVAR z-score
n_label   <- 12                 # labelled motifs per direction

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
register(SerialParam())         # multicore is not supported on Windows

lab_hiv  <- paste("Higher in", hiv_label)
lab_mock <- paste("Higher in", mock_label)
lab_ns   <- "Not significant"
cols <- c("#C43C39", "#3676A8", "grey75")
names(cols) <- c(lab_hiv, lab_mock, lab_ns)

# ---------------- Load object once ----------------
obj <- readRDS(obj_path)
DefaultAssay(obj) <- "ATAC"
md <- obj[[]]

unique_clusters <- sort(unique(as.character(md[[cluster_col]])))
message("Found clusters: ", paste(unique_clusters, collapse = ", "))

# ---------------- Motif matching once, on all peaks ----------------
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
)
id2name <- setNames(TFBSTools::name(pfm), names(pfm))

all_ranges <- granges(obj[["ATAC"]])
names(all_ranges) <- rownames(obj[["ATAC"]])
message("Matching motifs on ", length(all_ranges), " peaks (once)...")
motif_all <- matchMotifs(pfm, all_ranges, genome = BSgenome.Hsapiens.UCSC.hg38)

# ---------------- Per-cluster analysis ----------------
run_cluster <- function(cluster_id) {
  message("\n==================== Cluster ", cluster_id, " ====================")
  cdir <- file.path(out_dir, paste0("Cluster_", cluster_id))
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  
  cells_keep <- rownames(md)[
    as.character(md[[cluster_col]]) == cluster_id &
      as.character(md[[condition_col]]) %in% c(hiv_label, mock_label)
  ]
  cond_all <- table(as.character(md[cells_keep, condition_col]))
  print(cond_all)
  if (length(cond_all) < 2 || any(cond_all < min_cells_per_condition)) {
    message("Skipping cluster ", cluster_id, ": too few cells in one or both conditions.")
    return(NULL)
  }
  
  # counts for this cluster only
  counts_mat <- GetAssayData(obj, assay = "ATAC", layer = "counts")[, cells_keep, drop = FALSE]
  counts_mat <- counts_mat[, Matrix::colSums(counts_mat) > 0, drop = FALSE]
  counts_mat <- counts_mat[Matrix::rowSums(counts_mat > 0) >= min_cells_per_peak, , drop = FALSE]
  message("Using ", nrow(counts_mat), " peaks and ", ncol(counts_mat), " cells")
  
  idx <- match(rownames(counts_mat), names(all_ranges))
  se <- SummarizedExperiment(
    assays    = list(counts = counts_mat),
    rowRanges = all_ranges[idx],
    colData   = md[colnames(counts_mat), , drop = FALSE]
  )
  motif_ix <- motif_all[idx, ]
  
  # chromVAR
  se  <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)
  bg  <- getBackgroundPeaks(se)
  dev <- computeDeviations(object = se, annotations = motif_ix, background_peaks = bg)
  
  dev_scores <- deviationScores(dev)
  dev_scores[is.na(dev_scores)] <- 0
  ids <- rownames(dev_scores)
  rownames(dev_scores) <- make.unique(ifelse(ids %in% names(id2name), id2name[ids], ids))
  
  # HIV vs ssDNA
  cond    <- as.character(colData(se)[[condition_col]])
  cells_h <- colnames(dev_scores)[cond == hiv_label]
  cells_m <- colnames(dev_scores)[cond == mock_label]
  if (length(cells_h) < min_cells_per_condition || length(cells_m) < min_cells_per_condition) {
    message("Skipping cluster ", cluster_id, ": too few cells after filtering.")
    return(NULL)
  }
  x <- dev_scores[, cells_h, drop = FALSE]
  y <- dev_scores[, cells_m, drop = FALSE]
  
  p <- vapply(seq_len(nrow(x)), function(i) {
    tryCatch(wilcox.test(x[i, ], y[i, ])$p.value, error = function(e) NA_real_)
  }, numeric(1))
  
  res <- data.frame(
    cluster    = cluster_id,
    motif      = rownames(x),
    mean_HIV   = rowMeans(x),
    mean_ssDNA = rowMeans(y),
    delta      = rowMeans(x) - rowMeans(y),
    p_val      = p,
    p_adj      = p.adjust(p, method = "BH"),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      significance = case_when(
        p_adj < padj_cut & delta >=  min_delta ~ lab_hiv,
        p_adj < padj_cut & delta <= -min_delta ~ lab_mock,
        TRUE ~ lab_ns
      ),
      significance = factor(significance, levels = c(lab_hiv, lab_mock, lab_ns))
    ) %>%
    arrange(desc(delta))
  
  write.csv(res, file.path(cdir, paste0("cluster", cluster_id, "_chromVAR_HIV_vs_Mock.csv")),
            row.names = FALSE)
  message("Cells: ", length(cells_h), " HIV, ", length(cells_m), " Mock")
  print(table(res$significance))
  
  # volcano
  lab <- res %>%
    filter(significance != lab_ns) %>%
    group_by(significance) %>%
    slice_max(order_by = abs(delta), n = n_label, with_ties = FALSE) %>%
    ungroup()
  
  p_vol <- ggplot(res, aes(delta, -log10(pmax(p_adj, 1e-300)), colour = significance)) +
    geom_point(size = 1.2, alpha = 0.8) +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = c(-min_delta, min_delta), linetype = "dotted", linewidth = 0.4) +
    geom_text_repel(data = lab, aes(label = motif), size = 3, max.overlaps = Inf,
                    box.padding = 0.35, show.legend = FALSE) +
    scale_colour_manual(values = cols, drop = FALSE) +
    labs(
      title = paste0("Cluster ", cluster_id, ": motif activity, ", hiv_label, " vs ", mock_label),
      subtitle = paste0("chromVAR computed on cluster ", cluster_id, " cells only"),
      x = paste0("Difference in mean chromVAR z-score (", hiv_label, " - ", mock_label, ")"),
      y = expression(-log[10]("adj. P")), colour = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
  
  ggsave(file.path(cdir, paste0("cluster", cluster_id, "_volcano.png")),
         p_vol, width = 6.5, height = 5.5, dpi = 300)
  ggsave(file.path(cdir, paste0("cluster", cluster_id, "_volcano.pdf")),
         p_vol, width = 6.5, height = 5.5)
  
  # bar plot of top motifs
  bars <- res %>%
    filter(significance != lab_ns) %>%
    group_by(significance) %>%
    slice_max(order_by = abs(delta), n = n_label, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(motif = factor(motif, levels = motif[order(delta)]))
  
  if (nrow(bars) > 0) {
    p_bar <- ggplot(bars, aes(delta, motif, fill = significance)) +
      geom_col(width = 0.75) +
      geom_vline(xintercept = 0, linewidth = 0.3) +
      scale_fill_manual(values = cols, name = NULL, drop = FALSE) +
      labs(
        title = paste0("Cluster ", cluster_id, ": top differential motifs"),
        x = paste0("Difference in mean chromVAR z-score (", hiv_label, " - ", mock_label, ")"),
        y = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
    
    h <- max(4, 0.25 * nrow(bars) + 1.8)
    ggsave(file.path(cdir, paste0("cluster", cluster_id, "_top_motifs_bar.png")),
           p_bar, width = 6, height = h, dpi = 300)
    ggsave(file.path(cdir, paste0("cluster", cluster_id, "_top_motifs_bar.pdf")),
           p_bar, width = 6, height = h)
  }
  
  message("Finished cluster ", cluster_id)
  res
}

all_res <- list()
for (cl in unique_clusters) {
  r <- tryCatch(run_cluster(cl), error = function(e) {
    message("Cluster ", cl, " failed: ", conditionMessage(e))
    NULL
  })
  if (!is.null(r)) all_res[[cl]] <- r
  gc()
}

# ---------------- Combined results across clusters ----------------
if (length(all_res) > 0) {
  combined <- bind_rows(all_res)
  write.csv(combined, file.path(out_dir, "all_clusters_chromVAR_HIV_vs_Mock.csv"), row.names = FALSE)
  write.csv(
    combined %>% filter(significance != lab_ns) %>% arrange(cluster, desc(abs(delta))),
    file.path(out_dir, "all_clusters_chromVAR_significant.csv"), row.names = FALSE
  )
}

message("\nDone. Results saved in: ", out_dir)