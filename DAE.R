

# ============================================================
# WHOLE-POPULATION DIFFERENTIAL ATAC ACCESSIBILITY
# HIV vs Mock | ssDNA vs Mock | HIV vs ssDNA
# Donor-aware pseudobulk + nearest-gene annotation + volcanoes
# ============================================================

library(Seurat)
library(Signac)
library(Matrix)
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(GenomicRanges)

# ------------------------------------------------------------
# 1. OUTPUT DIRECTORY
# ------------------------------------------------------------

outdir <- "E:/DAE_results2/Whole_population"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. CHECK OBJECT AND ATAC COUNTS
# ------------------------------------------------------------

stopifnot("ATAC" %in% Assays(multiome))

DefaultAssay(multiome) <- "ATAC"

# Join split ATAC layers if using Seurat v5
#multiome <- JoinLayers(multiome, assay = "ATAC")

peak_counts <- GetAssayData(
  multiome,
  assay = "ATAC",
  layer = "counts"
)

meta <- multiome[[]][colnames(peak_counts), , drop = FALSE]

stopifnot(all(c("condition", "donor") %in% colnames(meta)))

meta$condition <- as.character(meta$condition)
meta$donor <- as.character(meta$donor)

cat("\nCells per condition and donor:\n")
print(table(meta$condition, meta$donor))

# ------------------------------------------------------------
# 3. AGGREGATE RAW COUNTS BY DONOR + CONDITION
# ------------------------------------------------------------

meta$pb_group <- paste(meta$donor, meta$condition, sep = "__")

# Require at least 20 cells in a donor-condition group
min_cells <- 20

group_sizes <- table(meta$pb_group)
valid_groups <- names(group_sizes[group_sizes >= min_cells])

pb_counts <- do.call(cbind, lapply(valid_groups, function(g) {
  cells <- rownames(meta)[meta$pb_group == g]
  Matrix::rowSums(peak_counts[, cells, drop = FALSE])
}))

rownames(pb_counts) <- rownames(peak_counts)
colnames(pb_counts) <- valid_groups

pb_meta <- data.frame(
  group = valid_groups,
  donor = sub("__.*$", "", valid_groups),
  condition = sub("^.*__", "", valid_groups),
  stringsAsFactors = FALSE
)
rownames(pb_meta) <- pb_meta$group

pb_meta$donor <- factor(pb_meta$donor)

cat("\nPseudobulk samples:\n")
print(pb_meta)

write.csv(
  pb_meta,
  file.path(outdir, "Pseudobulk_sample_metadata.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. DONOR-AWARE DIFFERENTIAL ACCESSIBILITY FUNCTION
# ------------------------------------------------------------

run_DAE <- function(counts, sample_meta, condition1, condition2) {
  
  # Keep only the two requested conditions
  keep <- sample_meta$condition %in% c(condition1, condition2)
  
  md <- sample_meta[keep, , drop = FALSE]
  counts_sub <- counts[, rownames(md), drop = FALSE]
  
  # Retain only donors represented in BOTH conditions
  donor_condition <- table(md$donor, md$condition)
  
  if (!all(c(condition1, condition2) %in% colnames(donor_condition))) {
    message("Skipping ", condition1, " vs ", condition2,
            ": one condition is missing.")
    return(NULL)
  }
  
  paired_donors <- rownames(donor_condition)[
    donor_condition[, condition1] > 0 &
      donor_condition[, condition2] > 0
  ]
  
  if (length(paired_donors) < 3) {
    message("Skipping ", condition1, " vs ", condition2,
            ": fewer than 3 paired donors.")
    return(NULL)
  }
  
  keep <- md$donor %in% paired_donors
  md <- md[keep, , drop = FALSE]
  counts_sub <- counts_sub[, rownames(md), drop = FALSE]
  
  md$donor <- droplevels(factor(md$donor))
  md$condition <- factor(
    md$condition,
    levels = c(condition2, condition1)
  )
  
  # Paired model: adjust for donor effects
  design <- model.matrix(~ donor + condition, data = md)
  
  y <- DGEList(counts = as.matrix(counts_sub))
  
  keep_peaks <- filterByExpr(y, design = design)
  y <- y[keep_peaks, , keep.lib.sizes = FALSE]
  
  if (nrow(y) == 0) return(NULL)
  
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)
  
  fit <- glmQLFit(y, design, robust = TRUE)
  
  coef_name <- paste0("condition", condition1)
  
  if (!coef_name %in% colnames(design)) {
    stop("Could not find coefficient: ", coef_name)
  }
  
  test <- glmQLFTest(fit, coef = coef_name)
  result <- topTags(test, n = Inf, sort.by = "PValue")$table
  
  result$peak <- rownames(result)
  result$comparison <- paste(condition1, "vs", condition2)
  
  # Positive logFC = more accessible in condition1
  result
}

# ------------------------------------------------------------
# 5. RUN ALL THREE COMPARISONS
# ------------------------------------------------------------

comparisons <- list(
  c("HIV", "Mock"),
    c("HIV", "ssDNA"),
  c("ssDNA", "Mock")
)

dae_list <- list()

for (cmp in comparisons) {
  
  condition1 <- cmp[1]
  condition2 <- cmp[2]
  
  cat("\nRunning:", condition1, "vs", condition2, "\n")
  
  res <- run_DAE(
    counts = pb_counts,
    sample_meta = pb_meta,
    condition1 = condition1,
    condition2 = condition2
  )
  
  if (!is.null(res)) {
    
    file_name <- paste0(
      "DAE_", condition1, "_vs_", condition2, ".csv"
    )
    
    write.csv(
      res,
      file.path(outdir, file_name),
      row.names = FALSE
    )
    
    dae_list[[paste(condition1, condition2, sep = "_vs_")]] <- res
  }
}

if (length(dae_list) == 0) {
  stop("No comparison produced results. Check donor-condition coverage.")
}

dae <- bind_rows(dae_list)

write.csv(
  dae,
  file.path(outdir, "DAE_all_comparisons.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. ANNOTATE PEAKS WITH NEAREST GENES
# ------------------------------------------------------------

# Confirm gene annotation is present
gene_annotation <- Annotation(multiome[["ATAC"]])

if (is.null(gene_annotation) || length(gene_annotation) == 0) {
  
  warning(
    "ATAC gene annotation is missing. Add the correct genome annotation ",
    "before attempting nearest-gene mapping."
  )
  
  dae_annotated <- dae
  
} else {
  
  # Convert peak IDs such as chr1-150561784-150562579
  # into genomic ranges
  peak_parts <- strsplit(dae$peak, "-", fixed = TRUE)
  
  valid_peak <- vapply(
    peak_parts,
    function(x) length(x) == 3 &&
      !is.na(suppressWarnings(as.numeric(x[2]))) &&
      !is.na(suppressWarnings(as.numeric(x[3]))),
    logical(1)
  )
  
  dae_valid <- dae[valid_peak, , drop = FALSE]
  peak_parts_valid <- peak_parts[valid_peak]
  
  peak_gr <- GRanges(
    seqnames = vapply(peak_parts_valid, `[`, character(1), 1),
    ranges = IRanges(
      start = as.integer(vapply(
        peak_parts_valid, `[`, character(1), 2
      )),
      end = as.integer(vapply(
        peak_parts_valid, `[`, character(1), 3
      ))
    )
  )
  
  names(peak_gr) <- dae_valid$peak
  
  # Find the closest annotated gene
  nearest_genes <- ClosestFeature(
    object = multiome,
    regions = peak_gr
  )
  
  # Find the returned region and gene-name columns
  region_col <- intersect(
    c("query_region", "peak", "region"),
    colnames(nearest_genes)
  )[1]
  
  gene_col <- intersect(
    c("gene_name", "gene", "gene_id"),
    colnames(nearest_genes)
  )[1]
  
  if (is.na(region_col) || is.na(gene_col)) {
    stop(
      "Could not identify region/gene columns in ClosestFeature output. ",
      "Inspect colnames(nearest_genes)."
    )
  }
  
  nearest_table <- data.frame(
    peak = as.character(nearest_genes[[region_col]]),
    nearest_gene = as.character(nearest_genes[[gene_col]])
  ) %>%
    distinct(peak, .keep_all = TRUE)
  
  dae_annotated <- dae %>%
    left_join(nearest_table, by = "peak")
}

write.csv(
  dae_annotated,
  file.path(outdir, "DAE_all_comparisons_annotated.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. PREPARE VOLCANO PLOT DATA
# ------------------------------------------------------------

dae_plot <- dae_annotated %>%
  filter(
    !is.na(logFC),
    !is.na(FDR),
    !is.na(comparison)
  ) %>%
  mutate(
    neglog10FDR = -log10(pmax(FDR, 1e-300)),
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "More accessible",
      FDR < 0.05 & logFC < 0 ~ "Less accessible",
      TRUE ~ "Not significant"
    )
  )

# Select up to 5 significant peaks per comparison for labeling
label_col <- if ("nearest_gene" %in% colnames(dae_plot)) {
  "nearest_gene"
} else {
  "peak"
}

label_data <- dae_plot %>%
  filter(FDR < 0.05)

if (label_col == "nearest_gene") {
  label_data <- label_data %>%
    filter(!is.na(nearest_gene), nearest_gene != "")
} else {
  label_data <- label_data %>%
    filter(!is.na(peak), peak != "")
}

label_data <- label_data %>%
  filter(
    !grepl(
      "^(MALAT1|MT-|RPS|RPL)",
      .data[[label_col]],
      ignore.case = TRUE
    )
  ) %>%
  group_by(comparison) %>%
  slice_min(order_by = FDR, n = 20, with_ties = FALSE) %>%
  ungroup()

# ------------------------------------------------------------
# 8. CREATE PUBLICATION-QUALITY VOLCANO PLOTS
# ------------------------------------------------------------

p_volcano <- ggplot(
  dae_plot,
  aes(x = logFC, y = neglog10FDR, color = significance)
) +
  geom_point(size = 0.65, alpha = 0.65) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_text_repel(
    data = label_data,
    aes(
      x = logFC,
      y = neglog10FDR,
      label = .data[[label_col]]
    ),
    inherit.aes = FALSE,
    size = 2.8,
    max.overlaps = Inf,
    box.padding = 0.35,
    show.legend = FALSE
  ) +
  facet_wrap(~comparison, scales = "free_y") +
  scale_color_manual(values = c(
    "More accessible" = "#C43C39",
    "Less accessible" = "#3676A8",
    "Not significant" = "grey75"
  )) +
  labs(
    title = "Differential chromatin accessibility",
    x = expression(log[2] * " fold change"),
    y = expression(-log[10]("FDR")),
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p_volcano)

# Save figures
ggsave(
  file.path(outdir, "Whole_population_volcano.pdf"),
  p_volcano,
  width = 12,
  height = 5,
  device = "pdf"
)

ggsave(
  file.path(outdir, "Whole_population_volcano.tiff"),
  p_volcano,
  width = 12,
  height = 5,
  dpi = 300,
  compression = "lzw"
)

# ------------------------------------------------------------
# 9. EXPORT SIGNIFICANT PEAKS
# ------------------------------------------------------------

significant_peaks <- dae_annotated %>%
  filter(!is.na(FDR), FDR < 0.05)

write.csv(
  significant_peaks,
  file.path(outdir, "DAE_significant_peaks_FDR05.csv"),
  row.names = FALSE
)

cat("\nAnalysis complete.\n")
cat("Results saved in:", outdir, "\n")
cat("Total tested peak-comparison rows:", nrow(dae), "\n")
cat("Significant peak-comparison rows:", nrow(significant_peaks), "\n")
