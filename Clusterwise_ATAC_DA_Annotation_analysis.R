
library(Seurat)
library(Signac)
library(ggplot2)
library(ggrepel)
library(dplyr)

# -------------------------------
# 1. Settings
# -------------------------------
outdir <- "E:/HIV_vs_Mock_ATAC_all_celltypes"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
obj<-multiome
DefaultAssay(obj) <- "ATAC"

celltype_col <- "Marie_celltype"
condition_col <- "condition"

min_cells <- 10
fdr_cut <- 0.05
lfc_cut <- 0.5
n_labels <- 20

# Check metadata and motif annotations
stopifnot(celltype_col %in% colnames(obj[[]]))
stopifnot(condition_col %in% colnames(obj[[]]))

print(table(obj[[condition_col]][, 1], useNA = "ifany"))
print(Assays(obj))

# Find all annotated cell types
cell_types <- sort(unique(na.omit(obj[[celltype_col]][, 1])))

# Confirm motif annotations exist before enrichment
if (is.null(Signac::Motifs(obj))) {
  stop("No motif annotations found. Add motif annotations before FindMotifs().")
}

# -------------------------------
# 2. Differential accessibility
#    HIV versus Mock per cell type
# -------------------------------
da_results <- list()

for (ct in cell_types) {
  
  message("Analysing: ", ct)
  
  cells <- colnames(obj)[
    !is.na(obj[[celltype_col]][, 1]) &
      obj[[celltype_col]][, 1] == ct &
      obj[[condition_col]][, 1] %in% c("HIV", "Mock")
  ]
  
  sub <- subset(obj, cells = cells)
  
  sub$condition <- factor(
    sub[[condition_col]][, 1],
    levels = c("Mock", "HIV")
  )
  
  counts <- table(sub$condition)
  print(counts)
  
  if (any(counts[c("Mock", "HIV")] < min_cells)) {
    message("Skipping ", ct, ": fewer than ", min_cells,
            " cells in one condition.")
    next
  }
  
  Idents(sub) <- "condition"
  
  da <- tryCatch(
    FindMarkers(
      sub,
      ident.1 = "HIV",
      ident.2 = "Mock",
      assay = "ATAC",
      test.use = "LR",
      latent.vars = "nCount_ATAC",
      min.pct = 0.05,
      logfc.threshold = 0
    ),
    error = function(e) {
      message("DA failed for ", ct, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(da) || nrow(da) == 0) next
  
  da$peak <- rownames(da)
  da$celltype <- ct
  
  da$direction <- ifelse(
    da$p_val_adj < fdr_cut & da$avg_log2FC > lfc_cut,
    "HIV increased",
    ifelse(
      da$p_val_adj < fdr_cut & da$avg_log2FC < -lfc_cut,
      "HIV decreased",
      "Not significant"
    )
  )
  
  da_results[[ct]] <- da
  
  # Save differential accessibility results
  safe_ct <- gsub("[^A-Za-z0-9]+", "_", ct)
  
  write.csv(
    da,
    file.path(outdir, paste0(safe_ct, "_DA_HIV_vs_Mock.csv")),
    row.names = FALSE
  )
  
  # Volcano plot
  da$neglog10FDR <- -log10(pmax(da$p_val_adj, 1e-300))
  
  label_df <- da |>
    filter(p_val_adj < fdr_cut) |>
    arrange(p_val_adj) |>
    head(n_labels)
  
  p <- ggplot(
    da,
    aes(x = avg_log2FC, y = neglog10FDR, color = direction)
  ) +
    geom_point(alpha = 0.7, size = 1.2) +
    geom_vline(
      xintercept = c(-lfc_cut, lfc_cut),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(fdr_cut),
      linetype = "dashed"
    ) +
    geom_text_repel(
      data = label_df,
      aes(label = peak),
      size = 2.5,
      max.overlaps = Inf
    ) +
    scale_color_manual(values = c(
      "HIV increased" = "#B2182B",
      "HIV decreased" = "#2166AC",
      "Not significant" = "grey75"
    )) +
    theme_classic(base_size = 12) +
    labs(
      title = paste("HIV vs Mock:", ct),
      x = "Log2 fold change in peak accessibility",
      y = "-Log10 adjusted P-value",
      color = NULL
    )
  
  ggsave(
    file.path(outdir, paste0(safe_ct, "_volcano.png")),
    p, width = 7, height = 5, dpi = 300
  )
}

# Combine all cell types
if (length(da_results) == 0) {
  stop("No DA results were produced. Check condition labels and cell counts.")
}

da_all <- bind_rows(da_results)

write.csv(
  da_all,
  file.path(outdir, "ALL_celltypes_DA_HIV_vs_Mock.csv"),
  row.names = FALSE
)

# -------------------------------
# 3. Motif enrichment
#    HIV-increased and decreased peaks
# -------------------------------
motif_results <- list()

for (ct in names(da_results)) {
  
  message("Motif enrichment: ", ct)
  
  da <- da_results[[ct]]
  safe_ct <- gsub("[^A-Za-z0-9]+", "_", ct)
  
  up_peaks <- da$peak[
    !is.na(da$p_val_adj) &
      da$p_val_adj < fdr_cut &
      da$avg_log2FC > lfc_cut
  ]
  
  down_peaks <- da$peak[
    !is.na(da$p_val_adj) &
      da$p_val_adj < fdr_cut &
      da$avg_log2FC < -lfc_cut
  ]
  
  for (direction in c("HIV_increased", "HIV_decreased")) {
    
    peaks <- if (direction == "HIV_increased") {
      up_peaks
    } else {
      down_peaks
    }
    
    if (length(peaks) < 5) {
      message("Skipping ", ct, " / ", direction,
              ": fewer than 5 significant peaks.")
      next
    }
    
    motifs <- tryCatch(
      FindMotifs(
        object = obj,
        features = peaks
      ),
      error = function(e) {
        message("Motif analysis failed for ", ct, " / ",
                direction, ": ", e$message)
        NULL
      }
    )
    
    if (is.null(motifs) || nrow(motifs) == 0) next
    
    motifs$celltype <- ct
    motifs$direction <- direction
    
    motif_results[[paste(ct, direction, sep = "__")]] <- motifs
    
    write.csv(
      motifs,
      file.path(
        outdir,
        paste0(safe_ct, "_", direction, "_motifs.csv")
      ),
      row.names = FALSE
    )
  }
}

if (length(motif_results) > 0) {
  motifs_all <- bind_rows(motif_results)
  
  write.csv(
    motifs_all,
    file.path(outdir, "ALL_celltypes_motif_enrichment.csv"),
    row.names = FALSE
  )
}

message("Analysis finished. Results saved in: ", outdir)
saveRDS(obj, "E:/HIV_vs_Mock_ATAC_all_celltypes/obj_HIV_Mock_ATAC.rds")

library(ggplot2)
library(ggrepel)
library(dplyr)
library(Signac)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86) # Adjust to your genome version if needed (e.g., hg38/GRCh38)

# --------------------------------------------------
# 1. Settings & Input Directories
# --------------------------------------------------
input_dir <- "E:/HIV_vs_Mock_ATAC_all_celltypes"
outdir <- file.path(input_dir, "annotated_cluster_plots")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load the combined master CSV file
master_csv_path <- file.path(input_dir, "ALL_celltypes_DA_HIV_vs_Mock.csv")
if (!file.exists(master_csv_path)) {
  stop("Master CSV file not found at: ", master_csv_path)
}

da_all <- read.csv(master_csv_path)

# Parameters for plotting
fdr_cut <- 0.05
lfc_cut <- 0.5
n_labels <- 10 # Number of top significant genes to label

# --------------------------------------------------
# 2. Annotate Peaks to Genes using Genomic Ranges
# --------------------------------------------------
message("Annotating peaks to nearest genes...")

# Convert peak string format (e.g., "chr11-6511000-6511500") to GenomicRanges
peak_coords <- unique(da_all$peak)
peak_gr <- Signac::StringToGRanges(peak_coords, sep = c("-", "-"))

# Extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC" # Ensure chromosome naming style matches (e.g.,chr11)

# Find the closest gene for each peak
closest_genes <- DistanceToNearest(peak_gr, annotations)
# Alternatively, use Signac's ClosestFeature if you have your Seurat object loaded:
# closest_genes <- ClosestFeature(obj, regions = peak_coords)

# Build a mapping data frame from peak to gene name
ann_df <- data.frame(
  peak = peak_coords,
  gene = annotations$gene_name[closest_genes$subjectHits],
  distance = closest_genes$distance
)

# Merge gene annotations back into the master data frame
da_all <- da_all %>% left_join(ann_df, by = "peak")

# --------------------------------------------------
# 3. Loop Through Each Cluster / Cell Type
# --------------------------------------------------
cell_types <- unique(da_all$celltype)

for (ct in cell_types) {
  message("Processing cluster: ", ct)
  
  # Subset data for the current cell type
  df_ct <- da_all %>% filter(celltype == ct)
  
  if (nrow(df_ct) == 0) next
  
  # Handle negative log10 p-value calculation safely
  df_ct$neglog10FDR <- -log10(pmax(df_ct$p_val_adj, 1e-300))
  
  # Select top significant genes to label on the volcano plot
  label_df <- df_ct %>%
    filter(p_val_adj < fdr_cut) %>%
    arrange(p_val_adj) %>%
    head(n_labels)
  
  # Generate Volcano Plot with Gene Labels
  p <- ggplot(df_ct, aes(x = avg_log2FC, y = neglog10FDR, color = direction)) +
    geom_point(alpha = 0.7, size = 1.2) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = label_df,
      aes(label = gene), # Labels points using the gene name!
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.5,
      force = 2
    ) +
    scale_color_manual(values = c(
      "HIV increased" = "#B2182B",
      "HIV decreased" = "#2166AC",
      "Not significant" = "grey75"
    )) +
    theme_classic(base_size = 12) +
    labs(
      title = paste("HIV vs Mock:", ct),
      subtitle = "Cluster-wise Differential Accessibility (Annotated)",
      x = "Log2 Fold Change (Accessibility)",
      y = "-Log10 Adjusted P-value",
      color = NULL
    )
  
  # Save individual cluster plot and table
  safe_ct <- gsub("[^A-Za-z0-9]+", "_", ct)
  
  ggsave(
    filename = file.path(outdir, paste0(safe_ct, "_annotated_volcano.png")),
    plot = p, width = 7, height = 5, dpi = 300
  )
  
  write.csv(
    df_ct,
    file.path(outdir, paste0(safe_ct, "_DA_annotated.csv")),
    row.names = FALSE
  )
}

# Save a complete master annotated table
write.csv(
  da_all,
  file.path(outdir, "ALL_celltypes_DA_HIV_vs_Mock_Annotated.csv"),
  row.names = FALSE
)

message("All cluster plots and annotated CSV files successfully saved in: ", outdir)
