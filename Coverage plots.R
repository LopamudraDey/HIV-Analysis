
library(Seurat)
library(Signac)
library(ggplot2)

multiome<-combined
# Check condition labels
table(multiome$condition)

# Coverage plot around IFI44L
CoveragePlot(
  object = multiome,
  region = "IFI44L",
  group.by = "condition",
  assay = "ATAC",
  annotation = TRUE,
  peaks = TRUE
)

library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
combined<-multiome
# Create a combined grouping variable
combined$cluster_condition <- paste(
  combined$seurat_clusters,
  combined$condition,
  sep = "_"
)

clusters <- sort(unique(as.character(combined$seurat_clusters)))

# Generate a coverage plot for each cluster
p_list <- lapply(clusters, function(cl) {
  
  obj_cl <- subset(
    combined,
    subset = seurat_clusters == cl
  )
  
  CoveragePlot(
    object = obj_cl,
    region = "IFI44L",
    group.by = "condition",
    assay = "ATAC",
    annotation = TRUE,
    peaks = TRUE
  ) +
    ggtitle(paste("Cluster", cl))
})

# Combine all cluster plots into one figure
all_clusters_plot <- wrap_plots(
  p_list,
  ncol = 2
)

all_clusters_plot

library(Signac)
library(Seurat)
library(ggplot2)

# Set the cluster identities
Idents(combined) <- "Marie_celltype"



p <- CoveragePlot(
  object = combined,
  region = "IFI44L",
  assay = "ATAC",
  group.by = "Marie_celltype",
  #split.by = "condition",
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 3000,
  extend.downstream = 3000
)
p
# Increase font sizes across all panels
p <- p & theme(
  text = element_text(size = 14),
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 14),
  plot.title = element_text(size = 16, face = "bold"),
  strip.text = element_text(size = 14, face = "bold")
)

p
# Save as a high-resolution PDF
ggsave(
  "IFI44L_ATAC_Coverage_All_Clusters.pdf",
  plot = all_clusters_plot,
  width = 14,
  height = max(6, 3 * ceiling(length(clusters) / 2))
)

genes <- c("IFI44L", "IFI6", "IFITM3", "IFI27")

plots <- lapply(genes, function(gene) {
  CoveragePlot(
    object = combined,
    region = gene,
    assay = "ATAC",
    group.by = "condition",
    annotation = TRUE,
    peaks = TRUE,
    extend.upstream = 3000,
    extend.downstream = 3000
  )
})
CoveragePlot(
  object = multiome,
  region = "IFI44",
  group.by = "condition",
  assay = "ATAC",
  annotation = TRUE,
  peaks = TRUE
)
names(plots) <- genes

# Display IFI44L
plots[["IFI44L"]]

###########################################HIV vs MOck##############################################################################
multiome<-readRDS("E:/multiome_with_Marie_celltype_atac.rds")
obj<-multiome
genes_to_plot <- c(
 "MX2","CCL4","GNLY","PLCG2","XAF1","TMSB10","NEAT1"
)
genes_to_plot <- c(
  "NEAT1"
)

genes_to_plot <- c(
  "TXN", "IFNG", "TMSB10"
)

genes_to_plot <- c(
  "PLCG2", "MX2", "TMSB10","IFI44L","XAF1","FAU","TRIM22","IFI44","Mx1", "ISG15"
)
genes_to_plot <- c(
  "PLCG2","CD2","IFI6","CD59","B2M","TMSB10"
)
genes_to_plot <- c(
  "CCL17","TMSB10","CCL3","CCL4","TXN","CCL22","GNLY"
)
genes_to_plot <- c(
  "IFI6","GZMB","CD52","CD3D","HLA-C","LAG3","IFNG","CD2"
)

genes_to_plot <- c(
  "IFI6","LAG3","TMA7","TMSB10"
)
genes_to_plot <- c(
  "PLCG2","IFI44L","IFI27","DDX5"
)
# 2. Metadata and condition labels

genes_to_plot <- c(
  "B2M", "CCL3", "CCL4", "PLCG2","IFI6"
)

genes_to_plot <- c(
  "B2M", "GZMB", "CCL3", "IFI6"
)

genes_to_plot <- c(
  "CCL17", "TXN"
)

genes_to_plot <- c(
  "NEAT1","MX2","RBMX","DDX5","NMT1","IGKC","IFITM1"
)
cluster_col <- "seurat_clusters"
condition_col <- "condition"

hiv_label <- "HIV"
mock_label <- "ssDNA"
cluster_id <- "9"

# 3. Output directory
output_dir <- "E:/CoveragePlots2/Cluster_9"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
obj<-multiome
# 4. Subset Cluster 0 and retain only HIV and Mock cells
metadata <- obj[[]]

cells_keep <- rownames(metadata)[
  as.character(metadata[[cluster_col]]) == cluster_id &
    as.character(metadata[[condition_col]]) %in%
    c(hiv_label, mock_label)
]

if (length(cells_keep) == 9) {
  stop("No HIV or Mock cells found in Cluster 6. Check metadata labels.")
}

sub_obj <- subset(obj, cells = cells_keep)

# Check cell numbers
print(table(sub_obj[[condition_col]][, 1]))

# 5. Check required assays
if (!"ATAC" %in% Assays(sub_obj)) {
  stop("ATAC assay not found. Check Assays(obj).")
}

if (!"RNA" %in% Assays(sub_obj)) {
  stop("RNA assay not found. Check Assays(obj).")
}

# 6. Generate coverage plots
DefaultAssay(sub_obj) <- "ATAC"

for (gene in genes_to_plot) {
  
  message("Plotting Cluster 2: ", gene)
  
  tryCatch({
    
    p <- CoveragePlot(
      object = sub_obj,
      region = gene,
      assay = "ATAC",
      features = gene,
      expression.assay = "RNA",
      group.by = condition_col,
      extend.upstream = 5000,
      extend.downstream = 5000,
      annotation = TRUE,
      peaks = TRUE,
      links = TRUE,
      show.bulk = FALSE
    )
    
    print(p)
    
    # Save PNG
    ggsave(
      filename = file.path(
        output_dir,
        paste0("Cluster_2_", gene, "_HIV_vs_ssDNA.png")
      ),
      plot = p,
      width = 12,
      height = 7,
      dpi = 300
    )
    
    # Save PDF
    ggsave(
      filename = file.path(
        output_dir,
        paste0("Cluster_2_", gene, "_HIV_vs_ssDNA.pdf")
      ),
      plot = p,
      width = 12,
      height = 7
    )
    
  }, error = function(e) {
    message(
      "Could not plot ", gene, ": ",
      conditionMessage(e)
    )
  })
}

message("Cluster 1 coverage plots completed.")
LinkPlot(
  object = multiome,
  region = "IFI44L",      # vector of gene symbols
  extend.upstream = 10000,
  extend.downstream = 10000
)
