
library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)

# 1. Read the CSV file
deg <- read_csv("E:/ssDNA_vs_Mock_GO/DE/Cluster_19_ssDNA_vs_Mock.csv")

# 2. Remove mitochondrial, Ribo, sex genes

sex_genes <- c("UTY", "USP9Y", "PRKY", "DDX3Y", "XIST", "TSIX", "EIF1AY","ZFY")

deg_filtered <- deg %>%
  filter(
    !grepl("^MT-", gene, ignore.case = TRUE),
    !grepl("^MTRNR", gene, ignore.case = TRUE),
    !grepl("MALAT1", gene, ignore.case = TRUE),
    !grepl("^RP[SL]", gene, ignore.case = TRUE),
    !gene %in% sex_genes  # Removes the specified sex-linked genes
  )
deg_filtered <- deg %>%
  filter(
    !grepl("^MT-", gene, ignore.case = TRUE),
    !grepl("^MTRNR", gene, ignore.case = TRUE),
    !grepl("MALAT1", gene, ignore.case = TRUE),
    !grepl("^RP[SL]", gene, ignore.case = TRUE)
  )
# 3. Define significance thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 0.25

# 4. Prepare data for volcano plot
deg_filtered <- deg_filtered %>%
  mutate(
    neg_log10_padj = -log10(pmax(p_val_adj, 1e-300)),
    regulation = case_when(
      p_val_adj < padj_cutoff & avg_log2FC > log2fc_cutoff ~ "Upregulated",
      p_val_adj < padj_cutoff & avg_log2FC < -log2fc_cutoff ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

# 5. Create volcano plot
volcano <- ggplot(
  deg_filtered,
  aes(x = avg_log2FC, y = neg_log10_padj, color = regulation)
) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not significant" = "grey70"
  )) +
  geom_vline(
    xintercept = c(-log2fc_cutoff, log2fc_cutoff),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype = "dashed"
  ) +
  labs(
    title = "Cluster 19: ssDNA vs Mock",
    subtitle = "",
    x = "Average log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Expression"
  ) +
  theme_classic(base_size = 13)

print(volcano)
top_genes <- deg_filtered %>%
  filter(
    p_val_adj < padj_cutoff,
    abs(avg_log2FC) > log2fc_cutoff
  ) %>%
  slice_min(p_val_adj, n = 20, with_ties = FALSE)

volcano <- volcano +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  )

volcano

library(Seurat)
library(ggplot2)

# 1. Ensure your active assay is set to RNA
DefaultAssay(obj) <- "RNA"


cluster_7_obj <- subset(obj, subset = seurat_clusters == "0") 

# 3. Generate the Seurat Violin Plot comparing ssDNA vs Mock
VlnPlot(cluster_7_obj, 
        features = c("PLCG2","TMSB10","IFI44L","IFITM1"),
        group.by = "condition",     # Change to your metadata column for treatment (e.g., "condition": ssDNA vs Mock)
        pt.size = 0,              # Size of individual cell points (set to 0 if too dense)
        combine = TRUE) 

# Switch to the ATAC assay to check chromatin accessibility peaks near these genes
DefaultAssay(obj) <- "ATAC"

# You can use CoveragePlot to visualize peaks around gene loci (e.g., B2M)
CoveragePlot(
  object = obj,
  region = "B2M", # Or the genomic coordinates of the peak
  group.by = "condition"
)

library(dplyr)
library(ggrepel)

# Select top 15 upregulated and top 15 downregulated genes
top_up <- deg_filtered %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.25
  ) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

top_down <- deg_filtered %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC < -0.25
  ) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

# Combine selected genes
top_genes <- bind_rows(top_up, top_down)

# Add labels to volcano plot
volcano <- volcano +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0
  )

print(volcano)

# Save the plot
ggsave(
  "Volcano_top15_up_down.png",
  plot = volcano,
  width = 10,
  height = 8,
  dpi = 300
)

# 6. Save the plot
ggsave(
  "Volcano_Cluster_0_HIV_vs_Mock_no_mito.png",
  plot = volcano,
  width = 8,
  height = 6,
  dpi = 300
)
