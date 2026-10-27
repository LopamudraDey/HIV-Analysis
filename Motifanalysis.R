library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

# ---------------------------------------------------------------
# 1. Load data and settings
# ---------------------------------------------------------------
obj <- readRDS("E:/multiome_with_Marie_celltype_atac.rds")

cluster_col   <- "seurat_clusters"
condition_col <- "condition"
cluster_id    <- "1"
hiv_label     <- "HIV"
mock_label    <- "ssDNA"

output_dir <- "E:/MotifEnrichment/Cluster_1"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------
# 2. Subset cluster 0, HIV + Mock cells only
# ---------------------------------------------------------------
md <- obj[[]]
cells_keep <- rownames(md)[
  as.character(md[[cluster_col]]) == cluster_id &
    as.character(md[[condition_col]]) %in% c(hiv_label, mock_label)
]
if (length(cells_keep) == 0) stop("No cells found - check cluster/condition labels.")

sub_obj <- subset(obj, cells = cells_keep)
print(table(sub_obj[[condition_col]][, 1]))

DefaultAssay(sub_obj) <- "ATAC"
Idents(sub_obj) <- condition_col

# ---------------------------------------------------------------
# 3. Add motif annotations (skip if your object already has them:
#    check with  Motifs(sub_obj)  )
# ---------------------------------------------------------------
if (is.null(Motifs(sub_obj))) {
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
  )
  # Chromosome names must match your peaks (UCSC style "chr1" for hg38)
  sub_obj <- AddMotifs(
    object  = sub_obj,
    genome  = BSgenome.Hsapiens.UCSC.hg38,
    pfm     = pfm
  )
}

# ---------------------------------------------------------------
# 4. Differentially accessible peaks: HIV vs Mock within cluster 0
# ---------------------------------------------------------------
da_peaks <- FindMarkers(
  object          = sub_obj,
  ident.1         = hiv_label,
  ident.2         = mock_label,
  test.use        = "LR",                 # recommended for ATAC
  latent.vars     = "nCount_ATAC",        # corrects for sequencing depth
  min.pct         = 0.05,
  logfc.threshold = 0,
  only.pos        = FALSE
)
da_peaks$peak <- rownames(da_peaks)
write.csv(da_peaks, file.path(output_dir, "DA_peaks_HIV_vs_Mock_cluster1.csv"), row.names = FALSE)

open_in_hiv  <- da_peaks %>% filter(p_val < 0.05, avg_log2FC >  0.25) %>% pull(peak)
open_in_mock <- da_peaks %>% filter(p_val < 0.05, avg_log2FC < -0.25) %>% pull(peak)
cat("Peaks more open in HIV :", length(open_in_hiv), "\n")
cat("Peaks more open in Mock:", length(open_in_mock), "\n")

# ---------------------------------------------------------------
# 5. Background peaks matched on GC content / accessibility
# ---------------------------------------------------------------
meta.feature <- GetAssayData(sub_obj, assay = "ATAC", layer = "meta.features")
# (use slot = "meta.features" on older Seurat/SeuratObject versions)

# Compute GC content etc. if missing
if (!"GC.percent" %in% colnames(meta.feature)) {
  sub_obj <- RegionStats(sub_obj, genome = BSgenome.Hsapiens.UCSC.hg38)
  meta.feature <- GetAssayData(sub_obj, assay = "ATAC", layer = "meta.features")
}

run_enrichment <- function(peaks_of_interest, label) {
  if (length(peaks_of_interest) < 10) {
    message("Too few peaks for ", label)
    return(NULL)
  }
  bg <- MatchRegionStats(
    meta.feature  = meta.feature,
    query.feature = meta.feature[peaks_of_interest, ],
    n             = 50000
  )
  enr <- FindMotifs(
    object   = sub_obj,
    features = peaks_of_interest,
    background = bg
  )
  enr$direction <- label
  write.csv(enr, file.path(output_dir, paste0("Motif_enrichment_", label, ".csv")), row.names = FALSE)
  enr
}

enr_hiv  <- run_enrichment(open_in_hiv,  "HIV_enriched")
enr_mock <- run_enrichment(open_in_mock, "Mock_enriched")

# ---------------------------------------------------------------
# 6. Plots
# ---------------------------------------------------------------
plot_top <- function(enr, label, n = 10) {
  if (is.null(enr)) return(invisible(NULL))
  top <- enr %>% filter(p.adjust < 0.05) %>% arrange(p.adjust) %>% head(n)
  if (nrow(top) == 0) { message("No significant motifs for ", label); return(invisible(NULL)) }
  
  p1 <- MotifPlot(sub_obj, motifs = top$motif, assay = "ATAC")
  ggsave(file.path(output_dir, paste0("Top_motifs_logo_", label, ".png")), p1, width = 10, height = 6)
  
  p2 <- ggplot(top, aes(x = reorder(motif.name, fold.enrichment), y = fold.enrichment, fill = -log10(p.adjust))) +
    geom_col() + coord_flip() + theme_classic() +
    labs(x = NULL, y = "Fold enrichment", title = paste("Cluster 0 -", label))
  ggsave(file.path(output_dir, paste0("Top_motifs_bar_", label, ".png")), p2, width = 7, height = 5)
}
plot_top(enr_hiv,  "HIV_enriched")
plot_top(enr_mock, "Mock_enriched")



# ============================================================
# Volcano plot from saved DA table: HIV vs Mock, cluster 0 (ATAC)
# Input : DA_peaks_HIV_vs_Mock_cluster0.csv  (from FindMarkers)
# ============================================================

library(dplyr)
library(ggplot2)
library(ggrepel)
library(Signac)
library(Seurat)

# ------------------------------------------------------------
# 0. SETTINGS
# ------------------------------------------------------------
output_dir  <- "E:/MotifEnrichment/Cluster_1"
da_file     <- file.path(output_dir, "DA_peaks_HIV_vs_Mock_cluster1.csv")

lfc_cut     <- 0.25     # same cutoffs you used for the motif input sets
p_cut       <- 0.05
use_adjusted <- FALSE   # FALSE: y = -log10(p_val) and colour by p_val (as in the motif input)
# TRUE : y = -log10(p_val_adj) and colour by p_val_adj
n_labels    <- 10       # labels per direction
max_dist_bp <- 50000    # only label peaks within this distance of a gene
annotate_genes <- TRUE  # needs the Seurat object (gene annotation) in memory

# ------------------------------------------------------------
# 1. READ TABLE
# ------------------------------------------------------------
da <- read.csv(da_file, stringsAsFactors = FALSE)
if (!"peak" %in% colnames(da)) da$peak <- da[[1]]

pcol <- if (use_adjusted) "p_val_adj" else "p_val"
da <- da %>%
  filter(!is.na(avg_log2FC), !is.na(.data[[pcol]])) %>%
  mutate(
    p_plot = .data[[pcol]],
    neglog10P = -log10(pmax(p_plot, 1e-300)),
    significance = case_when(
      p_plot < p_cut & avg_log2FC >  lfc_cut ~ "More accessible in HIV",
      p_plot < p_cut & avg_log2FC < -lfc_cut ~ "More accessible in Mock",
      TRUE ~ "Not significant"
    )
  )

cat("Peaks tested:", nrow(da), "\n")
print(table(da$significance))
cat("Peaks with adj. P <", p_cut, ":", sum(da$p_val_adj < p_cut, na.rm = TRUE), "\n")

# ------------------------------------------------------------
# 2. LABELS: nearest gene for the top peaks
# ------------------------------------------------------------
label_peaks <- da %>%
  filter(significance != "Not significant") %>%
  group_by(significance) %>%
  slice_min(order_by = p_plot, n = n_labels * 3, with_ties = FALSE) %>%   # extra, some get filtered
  ungroup()

label_peaks$label <- label_peaks$peak
if (annotate_genes) {
  if (!exists("obj")) obj <- readRDS("E:/multiome_with_Marie_celltype_atac.rds")
  DefaultAssay(obj) <- "ATAC"
  cf <- ClosestFeature(obj, regions = label_peaks$peak)
  region_col <- intersect(c("query_region", "peak", "region"), colnames(cf))[1]
  gene_col   <- intersect(c("gene_name", "gene", "gene_id"), colnames(cf))[1]
  nt <- data.frame(peak = as.character(cf[[region_col]]),
                   gene = as.character(cf[[gene_col]]),
                   distance = if ("distance" %in% colnames(cf)) cf$distance else NA_real_,
                   stringsAsFactors = FALSE) %>% distinct(peak, .keep_all = TRUE)
  label_peaks <- label_peaks %>% left_join(nt, by = "peak") %>%
    filter(!is.na(gene), gene != "", is.na(distance) | distance <= max_dist_bp) %>%
    mutate(label = gene)
}

label_peaks <- label_peaks %>%
  group_by(significance) %>%
  distinct(label, .keep_all = TRUE) %>%
  slice_min(order_by = p_plot, n = n_labels, with_ties = FALSE) %>%
  ungroup()

# ------------------------------------------------------------
# 3. VOLCANO
# ------------------------------------------------------------
n_hiv  <- sum(da$significance == "More accessible in HIV")
n_mock <- sum(da$significance == "More accessible in Mock")

p <- ggplot(da, aes(avg_log2FC, neglog10P, colour = significance)) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_hline(yintercept = -log10(p_cut), linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted", linewidth = 0.3) +
  geom_text_repel(data = label_peaks,
                  aes(x = avg_log2FC, y = neglog10P, label = label),
                  inherit.aes = FALSE, size = 3, max.overlaps = Inf,
                  box.padding = 0.4, segment.size = 0.3, show.legend = FALSE) +
  scale_colour_manual(values = c("More accessible in HIV"  = "#C43C39",
                                 "More accessible in Mock" = "#3676A8",
                                 "Not significant"         = "grey75")) +
  labs(title = "Cluster 0: HIV vs Mock (ATAC)",
       subtitle = paste0(n_hiv, " peaks more open in HIV, ", n_mock, " more open in Mock (",
                         pcol, " < ", p_cut, ", |log2FC| > ", lfc_cut, ")"),
       x = expression(log[2] * " fold change (HIV / Mock)"),
       y = bquote(-log[10] * "(" * .(if (use_adjusted) "adj. P" else "P") * ")"),
       colour = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

print(p)
ggsave(file.path(output_dir, "Volcano_cluster0_HIV_vs_Mock.png"), p, width = 6.5, height = 5.5, dpi = 300)

# ============================================================
# Motif enrichment plots from saved FindMotifs tables (cluster 0, HIV vs Mock)
# Inputs : Motif_enrichment_HIV_enriched.csv, Motif_enrichment_Mock_enriched.csv
# ============================================================

library(dplyr)
library(ggplot2)
library(ggrepel)

# ------------------------------------------------------------
# 0. SETTINGS
# ------------------------------------------------------------
output_dir   <- "E:/MotifEnrichment/Cluster_1"
padj_cut     <- 0.05
n_label      <- 10         # labelled motifs per direction in the volcano
n_dot        <- 10         # motifs per direction in the dot plot
drop_gc_rich <- FALSE      # TRUE hides KLF/SP/EGR/NRF1/ZBTB14-type GC-rich motifs
gc_rich_pat  <- "^(|WT1|PATZ1)"

# ------------------------------------------------------------
# 1. READ + COMBINE
# ------------------------------------------------------------
read_enr <- function(f, lab) {
  p <- file.path(output_dir, f)
  if (!file.exists(p)) { message("Missing: ", p); return(NULL) }
  d <- read.csv(p, stringsAsFactors = FALSE)
  d$direction <- lab
  d
}

enr <- bind_rows(
  read_enr("Motif_enrichment_HIV_enriched.csv",  "More open in HIV"),
  read_enr("Motif_enrichment_Mock_enriched.csv", "More open in Mock")
)
if (nrow(enr) == 0) stop("No motif tables found in ", output_dir)

enr <- enr %>%
  mutate(log2FE = log2(fold.enrichment),
         neglog10P = -log10(pmax(p.adjust, 1e-300)),
         significant = p.adjust < padj_cut)

if (drop_gc_rich) enr <- enr %>% filter(!grepl(gc_rich_pat, motif.name))

cat("Significant motifs (adj. P <", padj_cut, "):\n")
print(table(enr$direction, enr$significant))

write.csv(enr %>% filter(significant) %>% arrange(direction, p.adjust),
          file.path(output_dir, "Motif_enrichment_significant_combined.csv"), row.names = FALSE)

cols <- c("More open in HIV" = "#C43C39", "More open in Mock" = "#3676A8")

# ------------------------------------------------------------
# 2. VOLCANO-STYLE PLOT: enrichment vs significance
# ------------------------------------------------------------
lab <- enr %>% filter(significant) %>%
  group_by(direction) %>% slice_min(p.adjust, n = n_label, with_ties = FALSE) %>% ungroup()

p_vol <- ggplot(enr, aes(log2FE, neglog10P, colour = direction)) +
  geom_point(aes(alpha = significant), size = 1.6) +
  scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.25), guide = "none") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", linewidth = 0.4) +
  geom_text_repel(data = lab, aes(label = motif.name), size = 3,
                  max.overlaps = Inf, box.padding = 0.35, show.legend = FALSE) +
  scale_colour_manual(values = cols) +
  facet_wrap(~direction) +
  labs(title = "Cluster 1: motif enrichment in differentially accessible peaks",
       x = expression(log[2] * " fold enrichment (vs matched background)"),
       y = expression(-log[10]("adj. P")), colour = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none",
        strip.background = element_blank(), strip.text = element_text(face = "bold"))

print(p_vol)
ggsave(file.path(output_dir, "Motif_volcano_cluster1.pdf"), p_vol, width = 9, height = 5)
ggsave(file.path(output_dir, "Motif_volcano_cluster1.png"), p_vol, width = 9, height = 5, dpi = 300)


#############################################################Motif Enrichment Automate#############################################################################

library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

# ---------------------------------------------------------------
# 1. Global Settings & Data Loading
# ---------------------------------------------------------------
obj <- readRDS("E:/multiome_with_Marie_celltype_atac.rds")

cluster_col   <- "seurat_clusters"
condition_col <- "condition"
hiv_label     <- "HIV"
mock_label    <- "ssDNA"

base_output_dir <- "E:/MotifEnrichment2"
dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)

# Ensure default assay is ATAC
DefaultAssay(obj) <- "ATAC"

# Ensure motif annotations exist globally across the object once
if (is.null(Motifs(obj))) {
  message("Adding motif annotations to Seurat object...")
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
  )
  obj <- AddMotifs(
    object  = obj,
    genome  = BSgenome.Hsapiens.UCSC.hg38,
    pfm     = pfm
  )
}

# Find all unique clusters
md <- obj[[]]
unique_clusters <- sort(unique(as.character(md[[cluster_col]])))
message("Found clusters: ", paste(unique_clusters, collapse = ", "))

# ---------------------------------------------------------------
# 2. Main Loop Across All Clusters
# ---------------------------------------------------------------
for (cluster_id in unique_clusters) {
  message("\n========================================================")
  message("Processing Cluster: ", cluster_id)
  message("========================================================")
  
  output_dir <- file.path(base_output_dir, paste0("Cluster_", cluster_id))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Subset current cluster, HIV + Mock cells only
  cells_keep <- rownames(md)[
    as.character(md[[cluster_col]]) == cluster_id &
      as.character(md[[condition_col]]) %in% c(hiv_label, mock_label)
  ]
  
  if (length(cells_keep) < 20) {
    message("Skipping Cluster ", cluster_id, ": Too few cells (", length(cells_keep), ").")
    next
  }
  
  sub_obj <- subset(obj, cells = cells_keep)
  cond_counts <- table(sub_obj[[condition_col]][, 1])
  print(cond_counts)
  
  if (length(cond_counts) < 2 || any(cond_counts < 10)) {
    message("Skipping Cluster ", cluster_id, ": Insufficient cells in one or both conditions.")
    next
  }
  
  DefaultAssay(sub_obj) <- "ATAC"
  Idents(sub_obj) <- condition_col
  
  # ---------------------------------------------------------------
  # 3. Differentially Accessible Peaks: HIV vs Mock
  # ---------------------------------------------------------------
  da_peaks <- tryCatch({
    FindMarkers(
      object          = sub_obj,
      ident.1         = hiv_label,
      ident.2         = mock_label,
      test.use        = "LR",
      latent.vars     = "nCount_ATAC",
      min.pct         = 0.05,
      logfc.threshold = 0,
      only.pos        = FALSE
    )
  }, error = function(e) {
    message("FindMarkers failed for Cluster ", cluster_id, ": ", e$message)
    NULL
  })
  
  if (is.null(da_peaks) || nrow(da_peaks) == 0) {
    message("No DA peaks found for Cluster ", cluster_id, ". Skipping.")
    next
  }
  
  da_peaks$peak <- rownames(da_peaks)
  da_file <- file.path(output_dir, paste0("DA_peaks_HIV_vs_ssDNA_cluster_", cluster_id, ".csv"))
  write.csv(da_peaks, da_file, row.names = FALSE)
  
  lfc_cut <- 0.25
  p_cut <- 0.05
  
  open_in_hiv  <- da_peaks %>% filter(p_val < p_cut, avg_log2FC >  lfc_cut) %>% pull(peak)
  open_in_mock <- da_peaks %>% filter(p_val < p_cut, avg_log2FC < -lfc_cut) %>% pull(peak)
  cat("Peaks more open in HIV :", length(open_in_hiv), "\n")
  cat("Peaks more open in ssDNA:", length(open_in_mock), "\n")
  
  # ---------------------------------------------------------------
  # 4. Background Peaks & Motif Enrichment
  # ---------------------------------------------------------------
  meta.feature <- GetAssayData(sub_obj, assay = "ATAC", layer = "meta.features")
  if (!"GC.percent" %in% colnames(meta.feature)) {
    sub_obj <- RegionStats(sub_obj, genome = BSgenome.Hsapiens.UCSC.hg38)
    meta.feature <- GetAssayData(sub_obj, assay = "ATAC", layer = "meta.features")
  }
  
  run_enrichment <- function(peaks_of_interest, label) {
    if (length(peaks_of_interest) < 10) {
      message("Too few peaks for ", label)
      return(NULL)
    }
    bg <- MatchRegionStats(
      meta.feature  = meta.feature,
      query.feature = meta.feature[peaks_of_interest, ],
      n             = 50000
    )
    enr <- FindMotifs(
      object   = sub_obj,
      features = peaks_of_interest,
      background = bg
    )
    enr$direction <- label
    write.csv(enr, file.path(output_dir, paste0("Motif_enrichment_", label, ".csv")), row.names = FALSE)
    enr
  }
  
  enr_hiv  <- run_enrichment(open_in_hiv,  "HIV_enriched")
  enr_mock <- run_enrichment(open_in_mock, "ssDNA_enriched")
  
  # ---------------------------------------------------------------
  # 5. Volcano Plot (DA Peaks)
  # ---------------------------------------------------------------
  use_adjusted <- FALSE
  pcol <- if (use_adjusted) "p_val_adj" else "p_val"
  
  da_peaks <- da_peaks %>%
    filter(!is.na(avg_log2FC), !is.na(.data[[pcol]])) %>%
    mutate(
      p_plot = .data[[pcol]],
      neglog10P = -log10(pmax(p_plot, 1e-300)),
      significance = case_when(
        p_plot < p_cut & avg_log2FC >  lfc_cut ~ "More accessible in HIV",
        p_plot < p_cut & avg_log2FC < -lfc_cut ~ "More accessible in ssDNA",
        TRUE ~ "Not significant"
      )
    )
  
  # Nearest gene mapping for top labels
  n_labels <- 10
  label_peaks <- da_peaks %>%
    filter(significance != "Not significant") %>%
    group_by(significance) %>%
    slice_min(order_by = p_plot, n = n_labels * 3, with_ties = FALSE) %>%
    ungroup()
  
  label_peaks$label <- label_peaks$peak
  cf <- ClosestFeature(obj, regions = label_peaks$peak)
  region_col <- intersect(c("query_region", "peak", "region"), colnames(cf))[1]
  gene_col   <- intersect(c("gene_name", "gene", "gene_id"), colnames(cf))[1]
  
  if (!is.null(region_col) && !is.null(gene_col)) {
    nt <- data.frame(
      peak = as.character(cf[[region_col]]),
      gene = as.character(cf[[gene_col]]),
      distance = if ("distance" %in% colnames(cf)) cf$distance else NA_real_,
      stringsAsFactors = FALSE
    ) %>% distinct(peak, .keep_all = TRUE)
    
    label_peaks <- label_peaks %>% left_join(nt, by = "peak") %>%
      filter(!is.na(gene), gene != "", is.na(distance) | distance <= 50000) %>%
      mutate(label = gene)
  }
  
  label_peaks <- label_peaks %>%
    group_by(significance) %>%
    distinct(label, .keep_all = TRUE) %>%
    slice_min(order_by = p_plot, n = n_labels, with_ties = FALSE) %>%
    ungroup()
  
  n_hiv_count  <- sum(da_peaks$significance == "More accessible in HIV")
  n_mock_count <- sum(da_peaks$significance == "More accessible in ssDNA")
  
  p_volcano <- ggplot(da_peaks, aes(avg_log2FC, neglog10P, colour = significance)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted", linewidth = 0.3) +
    geom_text_repel(
      data = label_peaks,
      aes(x = avg_log2FC, y = neglog10P, label = label),
      inherit.aes = FALSE, size = 3, max.overlaps = Inf,
      box.padding = 0.4, segment.size = 0.3, show.legend = FALSE
    ) +
    scale_colour_manual(values = c(
      "More accessible in HIV"  = "#C43C39",
      "More accessible in Mock" = "#3676A8",
      "Not significant"         = "grey75"
    )) +
    labs(
      title = paste("Cluster", cluster_id, ": HIV vs ssDNA (ATAC)"),
      subtitle = paste0(n_hiv_count, " peaks more open in HIV, ", n_mock_count, " more open in ssDNA"),
      x = expression(log[2] * " fold change (HIV / ssDNA)"),
      y = bquote(-log[10] * "(" * .(if (use_adjusted) "adj. P" else "P") * ")"),
      colour = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
  
  ggsave(file.path(output_dir, paste0("Volcano_cluster_", cluster_id, "_HIV_vs_ssDNA.png")), p_volcano, width = 6.5, height = 5.5, dpi = 300)
  
  # ---------------------------------------------------------------
  # 6. Motif Enrichment Volcano Plots
  # ---------------------------------------------------------------
  read_enr <- function(filename, lab) {
    p <- file.path(output_dir, filename)
    if (!file.exists(p)) return(NULL)
    d <- read.csv(p, stringsAsFactors = FALSE)
    d$direction <- lab
    d
  }
  
  enr_combined <- bind_rows(
    read_enr("Motif_enrichment_HIV_enriched.csv",  "More open in HIV"),
    read_enr("Motif_enrichment_ssDNA_enriched.csv", "More open in ssDNA")
  )
  
  if (!is.null(enr_combined) && nrow(enr_combined) > 0) {
    padj_cut <- 0.05
    enr_combined <- enr_combined %>%
      mutate(
        log2FE = log2(fold.enrichment),
        neglog10P = -log10(pmax(p.adjust, 1e-300)),
        significant = p.adjust < padj_cut
      )
    
    write.csv(
      enr_combined %>% filter(significant) %>% arrange(direction, p.adjust),
      file.path(output_dir, "Motif_enrichment_significant_combined.csv"),
      row.names = FALSE
    )
    
    cols <- c("More open in HIV" = "#C43C39", "More open in Mock" = "#3676A8")
    lab_motifs <- enr_combined %>%
      filter(significant) %>%
      group_by(direction) %>%
      slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
      ungroup()
    
    p_motif_vol <- ggplot(enr_combined, aes(log2FE, neglog10P, colour = direction)) +
      geom_point(aes(alpha = significant), size = 1.6) +
      scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.25), guide = "none") +
      geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", linewidth = 0.4) +
      geom_text_repel(
        data = lab_motifs, aes(label = motif.name), size = 3,
        max.overlaps = Inf, box.padding = 0.35, show.legend = FALSE
      ) +
      scale_colour_manual(values = cols) +
      facet_wrap(~direction) +
      labs(
        title = paste("Cluster", cluster_id, ": Motif Enrichment"),
        x = expression(log[2] * " fold enrichment (vs matched background)"),
        y = expression(-log[10]("adj. P")), colour = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")
      )
    
    ggsave(file.path(output_dir, paste0("Motif_volcano_cluster_", cluster_id, ".pdf")), p_motif_vol, width = 9, height = 5)
    ggsave(file.path(output_dir, paste0("Motif_volcano_cluster_", cluster_id, ".png")), p_motif_vol, width = 9, height = 5, dpi = 300)
  }
  
  message("Finished processing Cluster ", cluster_id)
}

message("\nAll clusters processed successfully! Results saved under: ", base_output_dir)
  
  

# ------------------------------------------------------------
# 4. SUMMARY TABLES
# ------------------------------------------------------------
if (length(summary_list) > 0) {
  summary_df <- bind_rows(summary_list)
  write.csv(summary_df, file.path(output_dir, "Summary_all_clusters.csv"), row.names = FALSE)
  print(summary_df)
}
if (length(motif_sig_list) > 0)
  write.csv(bind_rows(motif_sig_list), file.path(output_dir, "Motif_enrichment_significant_all_clusters.csv"),
            row.names = FALSE)


cat("\nDone. Results in:", output_dir, "\n")