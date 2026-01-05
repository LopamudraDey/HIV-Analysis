.libPaths("D:/HIV3")
dir.create("D:/HIV3", showWarnings = FALSE)
install.packages("BiocManager")
.libPaths()

BiocManager::install(c(
  "BiocGenerics","IRanges","S4Vectors","GenomicRanges",
  "SummarizedExperiment","AnnotationDbi","Biobase",
  "ensembldb","GenomicFeatures","EnsDb.Hsapiens.v86",
  "BSgenome.Hsapiens.UCSC.hg38",
  "SeuratObject","Seurat","Signac"
), lib = "D:/HIV3", force = TRUE)

.libPaths("D:/HIV3")

library(GenomicRanges)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Signac)
library(Seurat)

install.packages("SoupX")
library(SoupX)
library(Seurat)
library(SoupX)
install.packages("hdf5r")
library(hdf5r)
library(hdf5r)
library(dplyr)
library(ggplot2) # For visualization
base <- "D:/sens2025621/ARCdata"
base
samples <- list(
  D2_HIV   = file.path(base, "D2_HIV_arc/outs"),
  D3_HIV   = file.path(base, "D3_HIV_arc/outs"),
  D2_Mock  = file.path(base, "D2_Mock_arc/outs"),
  D3_Mock  = file.path(base, "D3_Mock_arc/outs"),
  D2_ssDNA = file.path(base, "D2_ssDNA_arc/outs"),
  D3_ssDNA = file.path(base, "D3_ssDNA_arc/outs")
)
# --- SET YOUR FILE PATHS HERE ---
# Path to the folder containing 'matrix.mtx', 'barcodes.tsv', and 'features.tsv' for ALL droplets
D2_HIV_raw_data_dir <- "D:/sens2025621/ARCdata/D2_HIV_arc/outs/raw_feature_bc_matrix/" 
# Path to the folder containing 'matrix.mtx', 'barcodes.tsv', and 'features.tsv' for FILTERED cells
D2_HIV_filtered_data_dir <- "D:/sens2025621/ARCdata/D2_HIV_arc/outs/filtered_feature_bc_matrix/" 

# --------------------------------

# Load the raw (Table of Droplets) and filtered (Table of Counts) matrices
tod <- Read10X(data.dir = D2_HIV_raw_data_dir) # Table of Droplets (All UMI)
toc <- Read10X(data.dir = D2_HIV_filtered_data_dir) # Table of Counts (Filtered Cells)
tod_rna <- tod[["Gene Expression"]]
toc_rna <- toc[["Gene Expression"]]
# Check that the dimensions are compatible
print(paste("Raw Droplets (tod) dimensions:", paste(dim(tod), collapse = " x ")))
print(paste("Filtered Cells (toc) dimensions:", paste(dim(toc), collapse = " x ")))
seurat_obj <- CreateSeuratObject(counts = toc_rna, project = "D2_HIV_Ambient_Correction")

# Ensure they have the same feature/gene set, which they usually do for 10X output
# Set RNA as active assay
DefaultAssay(seurat_obj) <- "RNA"
# 1. Set the assay
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
# 2. Run the detection and save the output to a vector

# Find variable genes
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
# This will return a character vector of 2000 genes if successful
head(VariableFeatures(seurat_obj[["RNA"]]))
# Scale only the variable features
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

# Run PCA
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)

# Clustering and Dimensionality Reduction (like the UMAP you showed)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5) # Use the resolution you determined to be best
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Store the cluster labels
initial_clusters <- seurat_obj@meta.data$seurat_clusters
names(initial_clusters) <- colnames(seurat_obj)



# 1. Create the SoupChannel object
sc <- SoupChannel(tod = tod_rna, toc = toc_rna)

# 2. Add the preliminary cluster information
sc <- setClusters(sc, initial_clusters)

sc <- autoEstCont(sc) 

# The estimated contamination fraction (rho) for the entire dataset:
print(paste("Estimated global contamination fraction (rho):", sc$soupProfile$rho))

# 4. Adjust the counts to remove the ambient RNA
corrected_counts <- adjustCounts(sc)

# 5. Create a new Seurat object with the corrected data
seurat_corrected <- CreateSeuratObject(counts = corrected_counts, 
                                       meta.data = seurat_obj@meta.data)

# Carry over the dimensional reduction (UMAP) coordinates from the original object
# This ensures we are comparing the same cells and map structures.
seurat_corrected[["umap"]] <- seurat_obj[["umap"]]

# Carry over the cluster identities
seurat_corrected@meta.data$seurat_clusters <- seurat_obj@meta.data$seurat_clusters
FeaturePlot(seurat_corrected, features = "CD3D")
seurat_corrected[["umap"]] <- seurat_obj[["umap"]]

# 3. Transfer the cluster identities
seurat_corrected@meta.data$seurat_clusters <- seurat_obj@meta.data$seurat_clusters

# 4. CRITICAL STEP: Re-run Normalization on the corrected data
# This step computes the 'data' layer, which FeaturePlot requires.
seurat_corrected <- NormalizeData(seurat_corrected)

# 5. CRITICAL STEP: Re-run Scaling on the corrected data
# This step computes the 'scale.data' layer, which is also standard practice.
# Note: You need to set the list of variable features again before scaling.
seurat_corrected <- FindVariableFeatures(seurat_corrected, selection.method = "vst", nfeatures = 2000)
seurat_corrected <- ScaleData(seurat_corrected)
seurat_corrected <- RunPCA(seurat_corrected, features = VariableFeatures(seurat_corrected), npcs = 30)

# 6. Re-plot S100a8
# This should now work without the "layer 'data' is not found" error.
FeaturePlot(seurat_corrected, features = "Cd3D", label = TRUE)
seurat_corrected

FeaturePlot(seurat_obj, features = c("CD3D", "CD3E"))
FeaturePlot(seurat_corrected, features = c("CD3D", "CD3E"))

DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") + ggtitle("Before SoupX")
DimPlot(seurat_corrected, reduction = "umap", group.by = "seurat_clusters") + ggtitle("After SoupX")
library(dplyr)

# Before SoupX
avg_before <- AverageExpression(seurat_obj, features = "CD3D")$RNA
# After SoupX
avg_after <- AverageExpression(seurat_corrected, features = "CD3D")$RNA

# Compare
avg_before
avg_after
sc$soupProfile$rho


D2_HIV_clean<-seurat_corrected
D2_HIV_clean
saveRDS(D2_HIV_clean,"D2_HIV_clean.rds")
##########################################################################################################################################################
samples <- list(
  D2_HIV  = list(raw="D:/sens2025621/ARCdata/D2_HIV_arc/outs/raw_feature_bc_matrix/",
                 filtered="D:/sens2025621/ARCdata/D2_HIV_arc/outs/filtered_feature_bc_matrix/"),
  D3_HIV  = list(raw="D:/sens2025621/ARCdata/D3_HIV_arc/outs/raw_feature_bc_matrix/",
                 filtered="D:/sens2025621/ARCdata/D3_HIV_arc/outs/filtered_feature_bc_matrix/"),
  D2_Mock = list(raw="D:/sens2025621/ARCdata/D2_Mock_arc/outs/raw_feature_bc_matrix/",
                 filtered="D:/sens2025621/ARCdata/D2_Mock_arc/outs/filtered_feature_bc_matrix/"),
  D3_Mock = list(raw="D:/sens2025621/ARCdata/D3_Mock_arc/outs/raw_feature_bc_matrix/",
                 filtered="D:/sens2025621/ARCdata/D3_Mock_arc/outs/filtered_feature_bc_matrix/"),
  D2_ssDNA= list(raw="D:/sens2025621/ARCdata/D2_ssDNA_arc/outs/raw_feature_bc_matrix/",
                 filtered="D:/sens2025621/ARCdata/D2_ssDNA_arc/outs/filtered_feature_bc_matrix/"),
  D3_ssDNA= list(raw="D:/sens2025621/ARCdata/D3_ssDNA_arc/outs/raw_feature_bc_matrix/",
                 filtered="D:/sens2025621/ARCdata/D3_ssDNA_arc/outs/filtered_feature_bc_matrix/")
)
library(Seurat)
library(SoupX)

process_sample <- function(sample_name, paths){
  
  # 1. Load raw & filtered matrices
  tod <- Read10X(data.dir = paths$raw)[["Gene Expression"]]
  toc <- Read10X(data.dir = paths$filtered)[["Gene Expression"]]
  
  # 2. Create initial Seurat object (filtered counts)
  seurat_obj <- CreateSeuratObject(counts = toc, project = sample_name)
  
  # 3. Standard preprocessing for clustering
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = 30)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  
  # 4. Create SoupChannel
  sc <- SoupChannel(tod = tod, toc = toc)
  initial_clusters <- seurat_obj@meta.data$seurat_clusters
  names(initial_clusters) <- colnames(seurat_obj)
  sc <- setClusters(sc, initial_clusters)
  
  # 5. Estimate contamination
  sc <- autoEstCont(sc)
  cat("Estimated rho for", sample_name, "=", sc$soupProfile$rho, "\n")
  
  # 6. Correct counts
  corrected_counts <- adjustCounts(sc)
  
  # 7. Create corrected Seurat object
  seurat_corrected <- CreateSeuratObject(counts = corrected_counts, meta.data = seurat_obj@meta.data)
  
  # 8. Re-run RNA preprocessing on corrected counts
  DefaultAssay(seurat_corrected) <- "RNA"
  seurat_corrected <- NormalizeData(seurat_corrected)
  seurat_corrected <- FindVariableFeatures(seurat_corrected, selection.method = "vst", nfeatures = 2000)
  seurat_corrected <- ScaleData(seurat_corrected, features = VariableFeatures(seurat_corrected))
  seurat_corrected <- RunPCA(seurat_corrected, features = VariableFeatures(seurat_corrected), npcs = 30)
  
  # Optional: transfer UMAP from original object
  if("umap" %in% names(seurat_obj@reductions)){
    seurat_corrected[["umap"]] <- seurat_obj[["umap"]]
  }
  
  # Save corrected object
  saveRDS(seurat_corrected, paste0(sample_name, "_clean.rds"))
  
  return(seurat_corrected)
}
corrected_samples <- list()

for(s in names(samples)){
  cat("Processing sample:", s, "\n")
  corrected_samples[[s]] <- process_sample(s, samples[[s]])
}
# Merge all 6 samples
combined_rna <- merge(
  corrected_samples$D2_HIV,
  y = list(
    corrected_samples$D3_HIV,
    corrected_samples$D2_Mock,
    corrected_samples$D3_Mock,
    corrected_samples$D2_ssDNA,
    corrected_samples$D3_ssDNA
  ),
  add.cell.ids = names(corrected_samples)
)

# Check
combined_rna
Assays(combined_rna)
table(combined_rna$sample)
library(Signac)
library(Seurat)
library(GenomeInfoDb)

# Load ATAC from filtered h5
counts <- Read10X_h5("D:/sens2025621/ARCdata/D2_HIV_arc/outs/filtered_feature_bc_matrix.h5")

atac_counts <- counts$Peaks

# Start from SoupX-corrected RNA object
seu <- corrected_samples$D2_HIV

# Add ATAC
seu[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_counts[, colnames(seu)],  # match barcodes
  sep = c(":", "-"),
  genome = "hg38"
)

DefaultAssay(seu) <- "RNA"
library(Seurat)
library(Signac)

base <- "D:/sens2025621/ARCdata"

sample_names <- names(corrected_samples)

multiome_samples <- list()

for (s in sample_names) {
  
  message("Processing ", s)
  
  # ---- 1. Load filtered multiome counts ----
  h5_path <- file.path(base, paste0(s, "_arc/outs/filtered_feature_bc_matrix.h5"))
  counts <- Read10X_h5(h5_path)
  
  atac_counts <- counts$Peaks
  
  # ---- 2. Get SoupX-corrected RNA Seurat object ----
  seu <- corrected_samples[[s]]
  
  # ---- 3. Match barcodes (CRITICAL) ----
  common_cells <- intersect(colnames(seu), colnames(atac_counts))
  
  seu <- subset(seu, cells = common_cells)
  atac_counts <- atac_counts[, common_cells]
  
  # ---- 4. Add ATAC assay ----
  seu[["ATAC"]] <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = "hg38"
  )
  
  # ---- 5. Metadata ----
  seu$sample <- s
  
  # ---- 6. Store ----
  multiome_samples[[s]] <- seu
}
multiome_samples$D2_Mock
Assays(multiome_samples$D2_HIV)

combinedHIV <- merge(
  multiome_samples[[1]],
  y = multiome_samples[-1],
  add.cell.ids = names(multiome_samples)
)

combinedHIV
Assays(combinedHIV)
table(combinedHIV$sample)
add_atac <- function(seu, sample_name, base_path) {
  
  counts <- Read10X_h5(file.path(base_path, sample_name, "outs", "filtered_feature_bc_matrix.h5"))
  atac_counts <- counts$Peaks
  
  frag_path <- file.path(base_path, sample_name, "outs", "atac_fragments.tsv.gz")
  
  seu[["ATAC"]] <- CreateChromatinAssay(
    counts = atac_counts[, colnames(seu)],
    sep = c(":", "-"),
    genome = "hg38",
    fragments = frag_path
  )
  
  return(seu)
}
base <- "D:/sens2025621/ARCdata"

corrected_samples$D2_HIV   <- add_atac(corrected_samples$D2_HIV,   "D2_HIV_arc",   base)
corrected_samples$D3_HIV   <- add_atac(corrected_samples$D3_HIV,   "D3_HIV_arc",   base)
corrected_samples$D2_Mock  <- add_atac(corrected_samples$D2_Mock,  "D2_Mock_arc",  base)
corrected_samples$D3_Mock  <- add_atac(corrected_samples$D3_Mock,  "D3_Mock_arc",  base)
corrected_samples$D2_ssDNA <- add_atac(corrected_samples$D2_ssDNA, "D2_ssDNA_arc", base)
corrected_samples$D3_ssDNA <- add_atac(corrected_samples$D3_ssDNA, "D3_ssDNA_arc", base)

combinedHIV <- merge(
  corrected_samples$D2_HIV,
  y = list(
    corrected_samples$D3_HIV,
    corrected_samples$D2_Mock,
    corrected_samples$D3_Mock,
    corrected_samples$D2_ssDNA,
    corrected_samples$D3_ssDNA
  ),
  add.cell.ids = names(corrected_samples)
)

annotations <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
######################################################
# 2. QC
#####################################################
# Step 2: Add RNA QC metrics
combinedHIV[["percent.mt"]] <- PercentageFeatureSet(combinedHIV, pattern = "^MT-")

# Step 3: Ensure ATAC assay is active
DefaultAssay(combinedHIV) <- "ATAC"
Annotation(combinedHIV) <- annotations
# Step 4: Compute nucleosome signal and TSS enrichment
combinedHIV <- NucleosomeSignal(combinedHIV)
combinedHIV <- TSSEnrichment(combinedHIV)

# Step 5: Add additional QC metrics for plotting
# total RNA counts
combinedHIV$nCount_RNA <- colSums(LayerData(combinedHIV[["RNA"]], layer = "counts"))
# total ATAC counts
combinedHIV$nCount_ATAC <- colSums(LayerData(combinedHIV[["ATAC"]], layer = "counts"))

# Step 6: Violin plots for QC
VlnPlot(
  combinedHIV,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"),
  ncol = 5,
  pt.size = 0,
  group.by = "orig.ident"
)

summary(combinedHIV$nCount_RNA)
saveRDS(combinedHIV,"combinedHIV_beforeQC.rds")
# Step 7 (optional): Filter low-quality cells
combinedHIV <- subset(
  combinedHIV,
  subset = nCount_RNA > 500 & nCount_RNA < 50000 &
    nCount_ATAC > 1000 & nCount_ATAC < 100000 &
    TSS.enrichment > 2 &
    nucleosome_signal < 4 &
    percent.mt < 20
)

ncol(combinedHIV)

saveRDS(combinedHIV, "HIV_ATAC_RNA_combined_after_QC.rds")
#combined1<-readRDS("D:/OneDrive - Linköpings universitet/Documents/HIV_ATAC_RNA_combined_QC.rds")

length(colnames(combinedHIV))       
Assays(combinedHIV)

nrow(combinedHIV[["RNA"]])   # number of RNA genes
nrow(combinedHIV[["ATAC"]])  # number of ATAC peaks
BiocManager::install("glmGamPoi")
library(glmGamPoi,lib.loc = "D:/HIV1")

DefaultAssay(combinedHIV) <- "RNA"
combinedHIV <- SCTransform(combinedHIV, verbose = TRUE)
combinedHIV <- RunPCA(combinedHIV, npcs = 30)
DefaultAssay(combinedHIV) <- "ATAC"
combinedHIV <- RunTFIDF(combinedHIV)
combinedHIV <- FindTopFeatures(combinedHIV)
combinedHIV <- RunSVD(combinedHIV)


combinedHIV <- FindMultiModalNeighbors(
  combinedHIV,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:30, 2:30)
)
Reductions(combinedHIV)
