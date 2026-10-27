## =============================================================================
## HIV / Mock / ssDNA Multiome Study - Full Pipeline
## Per-sample SoupX + ATAC -> Merge -> QC -> Harmony integration -> WNN
## =============================================================================
.libPaths()
library(Seurat)
library(Signac)
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
# 1. Install BiocManager from CRAN
install.packages("BiocManager")
install.packages("sparseMatrixStats")
install.packages('sparseMatrixStats')
BiocManager::install(c("GenomeInfoDb", "GenomicRanges", "IRanges", "S4Vectors", "Rsamtools"))
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

BiocManager::install("biovizBase")
install.packages("hdf5r")
install.packages(c("Signac"))
# 2. Run your BiocManager installation command
BiocManager::install(version = "3.14")
library(SoupX)
install.packages("SoupX")
install.packages("harmony")
library(harmony)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(Matrix)
install.packages("DoubletFinder")
library(DoubletFinder)
## =============================================================================
## 0. Build annotation ONCE, globally (fixes: `annotations` undefined bug)
## =============================================================================
annotation <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

## =============================================================================
## 1. Sample sheet
## =============================================================================
samples <- data.frame(
  sample = c("D4_HIV", "D4_Mock", "D4_ssDNA","D2_HIV", "D2_Mock", "D2_ssDNA", "D3_HIV", "D3_Mock", "D3_ssDNA"),
  path = c(
    "E:/correctedARC_data/D4_HIV_arc_corr",
    "E:/correctedARC_data/D4_Mock_arc_corr",
    "E:/correctedARC_data/D4_ssDNA_arc_corr",
    "E:/Cellranger_data/sens2025621/ARCdata/D2_HIV_arc",
    "E:/Cellranger_data/sens2025621/ARCdata/D2_Mock_arc",
    "E:/Cellranger_data/sens2025621/ARCdata/D2_ssDNA_arc",
    "E:/Cellranger_data/sens2025621/ARCdata/D3_HIV_arc",
    "E:/Cellranger_data/sens2025621/ARCdata/D3_Mock_arc",
    "E:/Cellranger_data/sens2025621/ARCdata/D3_ssDNA_arc"
  ),
  condition = c("HIV", "Mock", "ssDNA","HIV", "Mock", "ssDNA","HIV", "Mock", "ssDNA"),
  donor     = c("D4", "D4", "D4","D2", "D2", "D2","D3", "D3", "D3")   #  this is donor ID, not timepoint
)

samples

## =============================================================================
## 2. Per-sample processing function
## =============================================================================
process_multiome <- function(sample_path, sample_name, condition, donor) {
  
  cat("\nProcessing:", sample_name, "\n")
  
  ###############################
  # RNA loading
  ###############################
  raw_path      <- paste0(sample_path, "/outs/raw_feature_bc_matrix/")
  filtered_path <- paste0(sample_path, "/outs/filtered_feature_bc_matrix/")
  
  tod <- Read10X(raw_path)
  toc <- Read10X(filtered_path)
  
  tod <- tod[["Gene Expression"]]
  toc <- toc[["Gene Expression"]]
  
  common_genes <- intersect(rownames(tod), rownames(toc))
  tod <- tod[common_genes, ]
  toc <- toc[common_genes, ]
  
  ###############################
  # Create RNA object + per-sample QC (fix: QC now runs BEFORE SoupX clustering)
  ###############################
  obj <- CreateSeuratObject(toc, project = sample_name)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  obj$sample    <- sample_name
  obj$condition <- condition
  obj$donor     <- donor
  
  obj <- subset(obj, subset = nFeature_RNA > 300 & nCount_RNA > 500 & percent.mt < 20)
  cat("  Cells after pre-SoupX QC:", ncol(obj), "\n")
  
  ###############################
  # Pre-SoupX clustering (needed by SoupX)
  ###############################
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  
  ###############################
  # SoupX
  ###############################
  sc <- SoupChannel(tod = tod, toc = toc[, colnames(obj)])
  sc <- setClusters(sc, setNames(as.character(obj$seurat_clusters), colnames(obj)))
  sc <- autoEstCont(sc, forceAccept = TRUE)
  corrected <- adjustCounts(sc, roundToInt = TRUE)
  
  obj <- CreateSeuratObject(corrected, project = sample_name)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj$sample    <- sample_name
  obj$condition <- condition
  obj$donor     <- donor
  
  ###############################
  # ATAC
  ###############################
  arc_outs <- paste0(sample_path, "/outs/")
  h5file   <- paste0(arc_outs, "filtered_feature_bc_matrix.h5")
  
  counts      <- Read10X_h5(h5file)
  atac_counts <- counts[["Peaks"]]
  
  # Restrict to standard chromosomes (fix: was missing)
  grange_counts   <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  standard_chroms <- seqnames(grange_counts) %in% standardChromosomes(grange_counts)
  atac_counts     <- atac_counts[as.vector(standard_chroms), ]
  
  frag_file <- paste0(arc_outs, "atac_fragments.tsv.gz")
  
  ###############################
  # Match RNA and ATAC barcodes
  ###############################
  shared <- intersect(colnames(obj), colnames(atac_counts))
  cat("  Shared RNA/ATAC barcodes:", length(shared), "\n")
  
  atac_counts <- atac_counts[, shared]
  obj <- subset(obj, cells = shared)
  
  chrom_assay <- CreateChromatinAssay(
    counts     = atac_counts,
    fragments  = frag_file,
    annotation = annotation,     # global annotation object, defined once above
    sep        = c(":", "-"),
    min.cells  = 10
  )
  
  obj[["ATAC"]] <- chrom_assay
  
  cat("  Final cells for", sample_name, ":", ncol(obj), "\n")
  
  return(obj)
}


multiome.list <- list()

for (i in 1:nrow(samples)) {
  multiome.list[[i]] <- process_multiome(
    sample_path = samples$path[i],
    sample_name = samples$sample[i],
    condition   = samples$condition[i],
    donor       = samples$donor[i]
  )
}

saveRDS(multiome.list, "E:/All_3_Donors_HIVMockssDNAv1.rds")
lapply(multiome.list, ncol)

## =============================================================================
## 4. Merge all 9 samples
## =============================================================================
# Step 1: Merge the first 3 samples
multiome_part1 <- merge(
  x = multiome.list[[1]],
  y = multiome.list[2:4],
  add.cell.ids = samples$sample[1:4]
)

# Step 2: Pre-prefix samples 5-9 individually (multiome_part1's barcodes
# are already unique from step 1, so they don't need re-prefixing)
multiome.list[5:9] <- lapply(5:9, function(i) {
  RenameCells(multiome.list[[i]], add.cell.id = samples$sample[i])
})

# Step 3: Merge without add.cell.ids - all barcodes already globally unique
multiome1 <- merge(
  x = multiome_part1,
  y = multiome.list[5:9]
)

multiome<-multiome1
# Clean up intermediate object and clear memory
rm(multiome_part1)
gc()
# 2D Table: Condition vs Donor
table(multiome1$condition, multiome1$donor)

# 3D Table: Sample vs Condition vs Donor
table(multiome1$sample, multiome1$condition, multiome1$donor)
saveRDS(multiome1, "E:/All_3_donors_HIVMockssDNAmerge_v1.rds")

## =============================================================================
## 5. Join layers (fix: required before LayerData/counts access works post-merge
##    in Seurat v5 - merge() keeps each sample as a separate counts layer)
## =============================================================================
DefaultAssay(multiome) <- "RNA"
multiome <- JoinLayers(multiome)

## =============================================================================
## 6. Post-merge QC
## =============================================================================
DefaultAssay(multiome) <- "RNA"
multiome[["percent.mt"]] <- PercentageFeatureSet(multiome, pattern = "^MT-")

DefaultAssay(multiome) <- "ATAC"
Annotation(multiome) <- annotation      # fix: use the global `annotation`, not undefined `annotations`

multiome <- NucleosomeSignal(multiome)
multiome <- TSSEnrichment(multiome)

# nCount_RNA / nFeature_RNA / nCount_ATAC already exist from per-sample object
# creation and are preserved through merge() - no need to recompute via LayerData.

VlnPlot(
  multiome,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"),
  group.by = "sample",
  ncol = 5,
  pt.size = 0
)

# fix: filtered result now reused for every downstream step (previously discarded)
multiome <- subset(
  multiome,
  subset = nFeature_RNA > 300 &
    nCount_RNA   > 500 &
    percent.mt   < 20 &
    nCount_ATAC  > 1000 &
    TSS.enrichment  > 1 &
    nucleosome_signal < 2
)

cat("Cells after merged QC:", ncol(multiome), "\n")
table(multiome$condition, multiome$donor)

saveRDS(multiome, "E:/All_3_Donors_HIVMockssDNA_QCfiltered_v1.rds")
#multiome<-readRDS("E:/AllHIVMockssDNA_QCfiltered.rds")
###############################################################################
# DoubletFinder - split samples
###############################################################################

sample.list <- SplitObject(
  multiome,
  split.by = "orig.ident"
)

names(sample.list)

lapply(
  sample.list,
  ncol
)



################D4_HIV##########################
obj <- sample.list$D4_HIV

DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)

nExp_poi <- round(
  0.045 * ncol(obj)
)

obj <- doubletFinder(
  obj,
  PCs = 1:30,
  pN = 0.25,
  pK = 0.09,
  nExp = nExp_poi,
  sct = FALSE
)

df.col <- grep(
  "^DF.classifications",
  colnames(obj@meta.data),
  value = TRUE
)

table(obj@meta.data[[df.col]])

obj <- subset(
  obj,
  cells = rownames(
    obj@meta.data[
      obj@meta.data[[df.col]] == "Singlet",
      ,
      drop = FALSE
    ]
  )
)

sample.list$D4_HIV <- obj

#####D4_Mock###################

obj <- sample.list$D4_Mock

DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)

nExp_poi <- round(
  0.045 * ncol(obj)
)

obj <- doubletFinder(
  obj,
  PCs = 1:30,
  pN = 0.25,
  pK = 0.09,
  nExp = nExp_poi,
  sct = FALSE
)

df.col <- grep(
  "^DF.classifications",
  colnames(obj@meta.data),
  value = TRUE
)

table(obj@meta.data[[df.col]])

obj <- subset(
  obj,
  cells = rownames(
    obj@meta.data[
      obj@meta.data[[df.col]] == "Singlet",
      ,
      drop = FALSE
    ]
  )
)

sample.list$D4_Mock <- obj

#########D4_ssDNA###############

obj <- sample.list$D4_ssDNA

DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)

nExp_poi <- round(
  0.045 * ncol(obj)
)

obj <- doubletFinder(
  obj,
  PCs = 1:30,
  pN = 0.25,
  pK = 0.09,
  nExp = nExp_poi,
  sct = FALSE
)

df.col <- grep(
  "^DF.classifications",
  colnames(obj@meta.data),
  value = TRUE
)

table(obj@meta.data[[df.col]])

obj <- subset(
  obj,
  cells = rownames(
    obj@meta.data[
      obj@meta.data[[df.col]] == "Singlet",
      ,
      drop = FALSE
    ]
  )
)

sample.list$D4_ssDNA <- obj





obj <- sample.list$D2_HIV

DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)

# Expected doublets
nExp_poi <- round(
  0.045 * ncol(obj)
)

nExp_poi
obj <- doubletFinder(
  obj,
  PCs = 1:30,
  pN = 0.25,
  pK = 0.09,
  nExp = nExp_poi,
  sct = FALSE
)
grep(
  "^DF.classifications",
  colnames(obj@meta.data),
  value = TRUE
)
df.col <- grep(
  "^DF.classifications",
  colnames(obj@meta.data),
  value = TRUE
)

table(obj@meta.data[[df.col]])
obj <- subset(
  obj,
  cells = rownames(
    obj@meta.data[
      obj@meta.data[[df.col]] == "Singlet",
      ,
      drop = FALSE
    ]
  )
)

sample.list$D2_HIV <- obj

#D2_Mock
obj <- sample.list$D2_Mock

DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)

# Expected doublets
nExp_poi <- round(
  0.045 * ncol(obj)
)

nExp_poi
obj <- doubletFinder(
  obj,
  PCs = 1:30,
  pN = 0.25,
  pK = 0.09,
  nExp = nExp_poi,
  sct = FALSE
)
grep(
  "^DF.classifications",
  colnames(obj@meta.data),
  value = TRUE
)
df.col <- grep(
  "^DF.classifications",
  colnames(obj@meta.data),
  value = TRUE
)

table(obj@meta.data[[df.col]])
obj <- subset(
  obj,
  cells = rownames(
    obj@meta.data[
      obj@meta.data[[df.col]] == "Singlet",
      ,
      drop = FALSE
    ]
  )
)

sample.list$D2_Mock <- obj


# D2_ssDNA

obj <- sample.list$D2_ssDNA

DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)

# Expected doublets
nExp_poi <- round(
  0.045 * ncol(obj)
)

nExp_poi
obj <- doubletFinder(
  obj,
  PCs = 1:30,
  pN = 0.25,
  pK = 0.09,
  nExp = nExp_poi,
  sct = FALSE
)
grep(
  "^DF.classifications",
  colnames(obj@meta.data),
  value = TRUE
)
df.col <- grep(
  "^DF.classifications",
  colnames(obj@meta.data),
  value = TRUE
)

table(obj@meta.data[[df.col]])
obj <- subset(
  obj,
    obj@meta.data[
      obj@meta.data[[df.col]] == "Singlet",
      ,
      drop = FALSE
    ]
  )
)

sample.list$D2_ssDNA <- obj

# D3_HIV
obj <- sample.list$D3_HIV

DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)

# Expected doublets
nExp_poi <- round(
  0.045 * ncol(obj)
)

nExp_poi
obj <- doubletFinder(
  obj,
  PCs = 1:30,
  pN = 0.25,
  pK = 0.09,
  nExp = nExp_poi,
  sct = FALSE
)
grep(
  "^DF.classifications",
  colnames(obj@meta.data),
  value = TRUE
)
df.col <- grep(
  "^DF.classifications",
  colnames(obj@meta.data),
  value = TRUE
)

table(obj@meta.data[[df.col]])
obj <- subset(
  obj,
  cells = rownames(
    obj@meta.data[
      obj@meta.data[[df.col]] == "Singlet",
      ,
      drop = FALSE
    ]
  )
)

sample.list$D3_HIV <- obj


#D3_Mock
obj <- sample.list$D3_Mock

DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)

# Expected doublets
nExp_poi <- round(
  0.045 * ncol(obj)
)

nExp_poi
obj <- doubletFinder(
  obj,
  PCs = 1:30,
  pN = 0.25,
  pK = 0.09,
  nExp = nExp_poi,
  sct = FALSE
)
grep(
  "^DF.classifications",
  colnames(obj@meta.data),
  value = TRUE
)
df.col <- grep(
  "^DF.classifications",
  colnames(obj@meta.data),
  value = TRUE
)

table(obj@meta.data[[df.col]])
obj <- subset(
  obj,
  cells = rownames(
    obj@meta.data[
      obj@meta.data[[df.col]] == "Singlet",
      ,
      drop = FALSE
    ]
  )
)

sample.list$D3_Mock <- obj


#D3_ssDNA
obj <- sample.list$D3_ssDNA

DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)

# Expected doublets
nExp_poi <- round(
  0.045 * ncol(obj)
)

nExp_poi
obj <- doubletFinder(
  obj,
  PCs = 1:30,
  pN = 0.25,
  pK = 0.09,
  nExp = nExp_poi,
  sct = FALSE
)
grep(
  "^DF.classifications",
  colnames(obj@meta.data),
  value = TRUE
)
df.col <- grep(
  "^DF.classifications",
  colnames(obj@meta.data),
  value = TRUE
)

table(obj@meta.data[[df.col]])
obj <- subset(
  obj,
  cells = rownames(
    obj@meta.data[
      obj@meta.data[[df.col]] == "Singlet",
      ,
      drop = FALSE
    ]
  )
)

sample.list$D3_ssDNA <- obj

# Check

lapply(
  sample.list,
  ncol
)
# merge

multiome <- merge(
  sample.list[[1]],
  y = sample.list[2:9],
  add.cell.ids = names(sample.list)
)

table(multiome$orig.ident)

ncol(multiome)

table(
  multiome$condition,
  multiome$donor
)

saveRDS(
  multiome,
  "E:/All_3_Donors_HIVMockssDNA_DoubletFiltered_v1.rds"
)
library(future)
names(multiome@graphs)
plan(sequential)
options(future.globals.maxSize = 8 * 1024^3)
multiome<-readRDS("E:/All_3_Donors_HIVMockssDNA_DoubletFiltered_v1.rds")


getwd()
DefaultAssay(multiome) <- "RNA"
multiome <- JoinLayers(multiome)
## =============================================================================
## 7. RNA processing + Harmony integration (batch = sample)
## =============================================================================
DefaultAssay(multiome) <- "RNA"
multiome <- NormalizeData(multiome)
multiome <- FindVariableFeatures(multiome, nfeatures = 3000)
multiome <- ScaleData(multiome)
multiome <- RunPCA(multiome, npcs = 30, reduction.name = "pca")

DefaultAssay(multiome) <- "RNA"

multiome <- RunHarmony(
  multiome,
  group.by.vars = "donor",
  reduction = "pca",
  reduction.save = "harmony_rna"
)
DefaultAssay(multiome) <- "RNA"

multiome <- FindNeighbors(
  multiome,
  reduction = "harmony_rna",
  dims = 1:30,
  graph.name = c("rna_nn", "rna_snn"),
  verbose = TRUE
)


multiome <- FindClusters(
  multiome,
  graph.name = "rna_snn",
  resolution = 0.5
)
multiome <- RunUMAP(
  multiome,
  reduction = "harmony_rna",
  dims = 1:30,
  reduction.name = "umap.rna",
  reduction.key = "rnaUMAP_"
)
head(multiome@meta.data)
Reductions(multiome)
p_rna <- DimPlot(
  multiome,
  reduction = "umap.rna",
  group.by = "seurat_clusters",split.by = "condition",
  label = TRUE,
  repel = TRUE,
  label.size = 7,
  pt.size = 0.5
) +
  ggtitle("RNA-seq UMAP") +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  )

p_rna
VlnPlot(multiome)
colnames(multiome@meta.data)
table(multiome$donor)
DimPlot(
  multiome,
  reduction = "umap.rna",
  group.by = "donor"
)

## =============================================================================
## 8. ATAC processing + Harmony integration on LSI
## =============================================================================
DefaultAssay(multiome) <- "ATAC"

multiome <- RunTFIDF(multiome)
multiome <- FindTopFeatures(multiome, min.cutoff = "q0")
multiome <- RunSVD(multiome)


DefaultAssay(multiome) <- "ATAC"

multiome <- RunHarmony(
  multiome,
  group.by.vars = "donor",
  reduction.use = "lsi",
  reduction.save = "harmony_atac",
  project.dim   = FALSE
)

# Generate ATAC UMAP using Harmony-corrected LSI
multiome <- RunUMAP(
  multiome,
  reduction = "harmony_atac",
  dims = 2:30,
  reduction.name = "umap.atac",
  reduction.key = "atacUMAP_"
)

# Plot by condition
p_atac <- DimPlot(
  multiome,
  reduction = "umap.atac",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  label.size = 7,
  pt.size = 0.5
) +
  ggtitle("ATAC-seq UMAP") +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  )
p_atac
# Plot by donor
DimPlot(
  multiome,
  reduction = "umap.atac",
  group.by = "condition"
)
saveRDS(multiome, "E:/All_3_donors_HIVMockssDNA_RNA+ATAC.rds")
## =============================================================================
## 9. WNN joint clustering on Harmony-corrected embeddings
## =============================================================================
multiome <- FindMultiModalNeighbors(
  multiome,
  reduction.list = list("harmony_rna", "harmony_atac"),   # not "pca"/"lsi"
  dims.list = list(1:30, 2:30)
)

multiome <- RunUMAP(
  multiome,
  nn.name        = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key  = "wnnUMAP_"
)

multiome <- FindClusters(
  multiome,
  graph.name = "wsnn",
  algorithm  = 3,
  resolution = 0.4,
  verbose    = FALSE
)

table(multiome$donor)
table(multiome$sample)

DimPlot(
  multiome,
  reduction = "wnn.umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
)
#multiome <- FindClusters(multiome, graph.name = "wsnn", algorithm = 3, resolution = 0.3)
library(ggplot2)

DimPlot(
  multiome,
  reduction = "wnn.umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  label.size = 6
) +
  ggtitle("WNN UMAP: Multiome Cell Clusters") +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13)
  )
DimPlot(multiome, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE) + NoLegend()
DimPlot(multiome, reduction = "wnn.umap", group.by = "donor")
DimPlot(multiome, reduction = "wnn.umap", group.by = "sample")
cat("Joint WNN clusters:", nlevels(multiome$seurat_clusters), "\n")
DimPlot(multiome, reduction = "wnn.umap", group.by = "donor")
DimPlot(multiome, reduction = "wnn.umap", group.by = "condition")
DimPlot(multiome, reduction = "wnn.umap", group.by = "sample")
FeatureScatter(
  subset(multiome, condition == "HIV"),
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)
VlnPlot(multiome, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_ATAC", "TSS.enrichment"), group.by = "condition", pt.size = 0)
FeaturePlot(
  multiome,
  features = c("ISG15", "IFIT1", "MX1"),
  reduction = "wnn.umap",
  split.by = "condition"
)
FeaturePlot(subset(multiome, condition == "HIV"), features = c("ISG15", "IFIT1", "MX1"), reduction = "wnn.umap")
## =============================================================================
## 10. Visual sanity checks
## =============================================================================
DimPlot(multiome, reduction = "wnn.umap", label = TRUE) + ggtitle("WNN clusters")
DimPlot(multiome, reduction = "wnn.umap", group.by = "condition") + ggtitle("By condition")
DimPlot(multiome, reduction = "wnn.umap", group.by = "sample") + ggtitle("By sample")
DimPlot(multiome, reduction = "wnn.umap", group.by = "condition", split.by = "donor") + ggtitle("By donor (check batch mixing)")
DimPlot(multiome, reduction = "wnn.umap", group.by = "seurat_clusters", split.by = "donor") + ggtitle("By donor (check batch mixing)")

DimPlot(multiome, reduction = "wnn.umap", group.by = "seurat_clusters", split.by = "condition") + ggtitle("By Condition")
table(multiome$seurat_clusters, multiome$condition)
table(multiome$seurat_clusters, multiome$sample)

## =============================================================================
## 11. Save final integrated object - ready for downstream DE/DA/composition analysis
## =============================================================================
saveRDS(multiome, "E:/All_3_donors_HIVMockssDNA_WNN_final_v1.rds")

multiome<-readRDS("E:/All_3_Donors_HIVMockssDNA_WNN_final_singleR_v1.rds")
##########################################################################################################################################################

hiv<-readRDS("D:/OneDrive - Linköpings universitet/Documents/combinedhiv_celltype.rds")

hiv
DimPlot(multiome, reduction = "wnn.umap") + ggtitle("WNN clusters")
DimPlot(multiome, reduction = "wnn.umap", group.by = "condition") + ggtitle("By condition")
DimPlot(hiv, reduction = "wnnumap", group.by = "orig.ident") + ggtitle("By sample (check batch mixing)")
DimPlot(hiv, reduction = "wnnumap", group.by = "condition", split.by = "donor") + ggtitle("By donor")
DimPlot(hiv, reduction = "wnnumap", group.by = "condition", split.by = "orig.ident") 

multiome<-readRDS("E:/AllHIVMockssDNA_WNN_final.rds")
BiocManager::install(c("SingleR", "celldex", "SingleCellExperiment"))
library(SingleR)
library(celldex)
library(SingleCellExperiment)

DefaultAssay(multiome) <- "RNA"

multiome <- JoinLayers(
  multiome,
  assay = "RNA"
)
DefaultAssay(multiome) <- "RNA"

# check features
head(rownames(multiome))
Assays(multiome)
Layers(multiome[["RNA"]])
Layers(multiome[["ATAC"]])
sce <- as.SingleCellExperiment(multiome)

ref <- MonacoImmuneData()
pred <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main
)


multiome$SingleR <- pred$labels

DimPlot(
  multiome,
  reduction = "wnn.umap",
  group.by = "SingleR",
  split.by="condition",
  label = TRUE,
  repel = TRUE
)

DimPlot(
  multiome,
  reduction = "wnn.umap",
  group.by = "SingleR",
  split.by="donor",
  label = TRUE,
  repel = TRUE
)
DimPlot(
  multiome,
  reduction = "wnn.umap",
  group.by = "SingleR",
  label = TRUE,
  repel = TRUE
)
saveRDS(multiome, "E:/All_3_Donors_HIVMockssDNA_WNN_final_singleR_v1.rds")

######################################### Take 3 Donors #################################################
multiome<-readRDS("E:/All_4_Donors_HIVMockssDNA_WNN_final_singleR_v2.rds")
table(multiome$orig.ident)

# Subset to keep only non-D1 samples
multiome <- subset(multiome, subset = !grepl("^D1_", orig.ident))

# Verify the remaining identity levels and cell counts
table(multiome$orig.ident)

DimPlot(multiome, reduction = "wnn.umap", label = TRUE) + ggtitle("WNN clusters")
DimPlot(multiome, reduction = "wnn.umap", group.by = "condition") + ggtitle("By condition")
DimPlot(multiome, reduction = "wnn.umap", group.by = "sample") + ggtitle("By sample")
DimPlot(multiome, reduction = "wnn.umap", group.by = "condition", split.by = "donor") + ggtitle("By donor")

table(multiome$orig.ident, multiome$seurat_clusters)
##########################################################################################################
# 1. Re-calculate RNA variable features & PCA
multiome <- FindVariableFeatures(multiome, selection.method = "vst", nfeatures = 2000)
multiome <- ScaleData(multiome)
multiome <- RunPCA(multiome, npcs = 30)

# 2. Re-calculate ATAC LSI (if using Signac)
multiome <- FindTopFeatures(multiome, min.cutoff = 'q0')
multiome <- RunTFIDF(multiome)
multiome <- RunSVD(multiome)

# 3. Re-build WNN (Weighted Nearest Neighbor) graph & clusters
multiome <- FindMultiModalNeighbors(
  multiome, 
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:30, 2:30) # Skip LSI component 1
)
multiome <- RunUMAP(multiome, nn.name = "weighted.nn", reduction.name = "wnn.umap")
multiome <- FindClusters(multiome, graph.name = "wsnn", resolution = 0.3)


lineage_markers <- c(
  "CD3D", "CD3E", "CD4", "CD8A",      # T cells
  "NKG7", "GNLY", "KLRD1",             # NK
  "CD14", "LYZ", "FCGR3A",             # Monocytes
  "MS4A1", "CD79A", "CD19",            # B cells
  "PPBP", "PF4",                       # Platelets
  "HBB"                                # RBC/erythroid, if present
)
library(Seurat)
library(ggplot2)

## ============================================================
## T cell + DC marker dot plot across WNN clusters
## ============================================================
## Assumes `multiome` already has seurat_clusters from FindClusters()
## on the WNN graph (see step 11 of the pipeline script).

DefaultAssay(multiome) <- "RNA"

marker_groups <- list(
  "Pan T cell" = c("CD3D", "CD3E", "CD3G", "CD2", "TRAC", "TRBC2"),
  "CD4 T" = c("CD4", "IL7R", "CCR7"),
  "CD8 T" = c("CD8A", "CD8B", "GZMK"),
  "Treg" = c("FOXP3", "IL2RA", "CTLA4"),
  "Naive / Memory" = c("SELL", "TCF7", "CCR7", "IL7R"),
  "Cytotoxic / Effector" = c("GZMB", "GZMA", "PRF1", "NKG7", "IFNG"),
  "Pan DC" = c("ITGAX", "HLA-DRA", "HLA-DRB1"),
  "cDC1" = c("CLEC9A", "XCR1", "BATF3", "IRF8"),
  "cDC2" = c("CD1C", "FCER1A", "CLEC10A", "SIRPA"),
  "pDC" = c("LILRA4", "CLEC4C", "IRF7", "IL3RA"),
  "B cells" = c("MS4A1", "CD79A", "CD19")
)

## ============================================================
## Combine + filter to genes actually present in the object
## ============================================================

markers_all <- unique(unlist(marker_groups, use.names = FALSE))
markers <- markers_all[markers_all %in% rownames(multiome)]
missing_markers <- setdiff(markers_all, markers)

cat("\nMissing markers (not found in this object):\n"); print(missing_markers)
cat("\nMarkers used in plot:\n"); print(markers)

## ============================================================
## Marker-group boundary positions (auto-calculated)
## ============================================================

group_lengths <- sapply(marker_groups, function(x) sum(x %in% markers))
group_ends    <- cumsum(group_lengths)
group_starts  <- c(1, head(group_ends, -1) + 1)
group_centers <- (group_starts + group_ends) / 2
group_labels  <- names(marker_groups)

## ============================================================
## Build DotPlot data
## ============================================================

dp <- DotPlot(
  multiome,
  features = markers,
  group.by = "seurat_clusters",
  scale = TRUE
)$data

dp$features.plot <- factor(dp$features.plot, levels = markers)

n_levels <- length(levels(factor(dp$id)))
label_y  <- n_levels + 0.55
ylim_max <- n_levels + 0.5

## ============================================================
## Final plot (geom_vline removed)
## ============================================================

p <- ggplot(
  dp,
  aes(x = features.plot, y = id, size = pct.exp, fill = avg.exp.scaled)
) +
  geom_point(shape = 21, colour = "black", stroke = 0.35) +
  scale_fill_gradient(
    low = "yellow", high = "red",
    name = "Scaled\nexpression"
  ) +
  scale_size(range = c(0.5, 7), name = "Percent\nexpressed") +
  labs(x = NULL, y = "Cluster", title = "T cell and DC markers by cluster") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9, colour = "black"),
    axis.text.y = element_text(size = 11, face = "bold", colour = "black"),
    axis.ticks  = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 9),
    plot.margin  = margin(t = 15, r = 20, b = 50, l = 20),
    plot.title   = element_text(face = "bold", hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(0.5, ylim_max), clip = "off")

for (i in seq_along(group_centers)) {
  p <- p + annotate(
    "text", x = group_centers[i], y = label_y,
    label = group_labels[i], fontface = "bold", size = 3.5
  )
}

p

#DotPlot(multiome, features = marker_groups) + RotatedAxis()
library(Seurat)

# Define cell type annotations based on cluster IDs
cell_type_annotations <- c(
  "0"  = "Pan T cell (Treg)",
  "1"  = "Cytotoxic CD8 T cell",
  "2"  = "B cell",
  "3"  = "Pan T cell",
  "4"  = "CD8 T cell",
  "5"  = "Pan T cell",
  "6"  = "Naive/Memory CD8 T cell",
  "7"  = "Pan T cell (Treg)",
  "8"  = "CD4 T cell (Treg)",
  "9"  = "Pan T cell",
  "10" = "cDC2",
  "11" = "CD4 T cell / cDC2",
  "12" = "Cytotoxic / Effector T cell",
  "13" = "Cytotoxic / Memory T cell",
  "14" = "Pan DC (cDC1/cDC2)",
  "15" = "Pan DC (cDC1)",
  "16" = "Pan DC / pDC",
  "17" = "pDC / cDC1",
  "18" = "CD4 T cell / pDC",
  "19"="CD4+ T cells / Myeloid Doublets"
)

# Rename cluster identities in your Seurat object
multiome <- RenameIdents(multiome, cell_type_annotations)

# Store annotations in metadata
multiome$cell_type <- Idents(multiome)

# Visualize annotated clusters on your reduction plot (e.g., UMAP)
DimPlot(multiome, reduction = "wnn.umap", label = TRUE, repel = TRUE)


DimPlot(multiome, reduction = "wnn.umap", split.by ="donor", label = TRUE, repel = TRUE)
library(Seurat)
library(dplyr)

DefaultAssay(multiome) <- "RNA"
multiome <- JoinLayers(multiome)

# 1. Reset Idents to numerical clusters if changed
Idents(multiome) <- "seurat_clusters"

# 2. Find top markers for all clusters
all_markers <- FindAllMarkers(
  multiome,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# 3. Extract top 10 genes per cluster
top20 <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

write.csv(all_markers,"E:/All_3_Donors_Markers_All_Clusters.csv")
write.csv(top20,"F:/Elements/All_3_donors_Markers_top20.csv")

##########################HIV vs Mock Marker DOT PLOT #####################
library(Seurat)
library(dplyr)
library(ggplot2)

DefaultAssay(HIV1) <- "RNA"
cell_types <- unique(na.omit(as.character(HIV1$Marie_celltype)))

markers_by_celltype <- lapply(cell_types, function(ct) {
  
  cells_ct <- colnames(HIV1)[
    !is.na(HIV1$Marie_celltype) &
      HIV1$Marie_celltype == ct
  ]
  
  obj_ct <- subset(HIV1, cells = cells_ct)
  
  # Count cells in each condition
  condition_counts <- table(obj_ct$condition)
  
  # Skip cell types with fewer than 3 cells in either group
  if (!all(c("HIV", "Mock") %in% names(condition_counts)) ||
      condition_counts["HIV"] < 3 ||
      condition_counts["Mock"] < 3) {
    message("Skipping ", ct, ": fewer than 3 cells in HIV or Mock")
    return(NULL)
  }
  
  Idents(obj_ct) <- "condition"
  
  result <- FindMarkers(
    obj_ct,
    ident.1 = "HIV",
    ident.2 = "Mock",
    assay = "RNA",
    test.use = "wilcox",
    min.pct = 0.1,
    logfc.threshold = 0.25
  )
  
  result$gene <- rownames(result)
  result$celltype <- ct
  
  return(result)
})

# Combine results, excluding skipped cell types
markers_by_celltype <- bind_rows(markers_by_celltype)

# Check results
dim(markers_by_celltype)
head(markers_by_celltype)
write.csv(
  markers_by_celltype,
  "E:/DEGs_HIV_vs_Mock_all_celltypes.csv",
  row.names = FALSE
)
 table(HIV1$Marie_celltype, HIV1$condition)

counts[
  rowSums(counts[, intersect(c("HIV", "Mock"), colnames(counts)), drop = FALSE] < 3) > 0,
  ,
  drop = FALSE
]
# Top 5 genes increased in HIV
library(dplyr)
library(Seurat)
library(ggplot2)

library(Seurat)
library(dplyr)
library(ggplot2)


# Keep only HIV and Mock (check the exact spelling with table(HIV1$condition))
HIV1 <- subset(HIV1, subset = condition %in% c("Mock", "HIV"))
HIV1$condition <- factor(HIV1$condition, levels = c("Mock", "HIV"))
# Mitochondrial and ribosomal genes (works for human and mouse names)
exclude <- "^MT-|^MTRNR|^RP[SL][0-9]|^RPLP|^SNHG|^LINC00278|^UTY$|^USP9Y$|^PRKY$|^DDX3Y$|^KDM5D$|^EIF1AY$|^RPS4Y1$|^XIST$|^TSIX$"

top5_HIV <- markers_by_celltype %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0,
         !grepl(exclude, gene, ignore.case = TRUE)) %>%
  group_by(celltype) %>%
  slice_max(order_by = avg_log2FC, n = 2, with_ties = FALSE) %>%
  ungroup()

top5_HIV <- markers_by_celltype %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5,
         pct.1 > 0.25, pct.2 < 0.3,
         !grepl(exclude, gene, ignore.case = TRUE)) %>%
  group_by(celltype) %>%
  slice_max(order_by = pct.1 - pct.2, n = 2, with_ties = FALSE) %>%
  ungroup()
# Check results
top5_HIV %>% select(celltype, gene, avg_log2FC, p_val_adj)

# Unique genes, kept in cell-type order (gives the diagonal pattern)
genes_to_plot <- unique(top5_HIV$gene)
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(HIV1)]

# ---- 2. Make sure cell names are unique (fixes duplicate row.names error) ----
if (anyDuplicated(colnames(HIV1)) > 0) {
  HIV1 <- RenameCells(HIV1, new.names = make.unique(colnames(HIV1)))
}

# ---- 3. Identities: cell type | condition -----------------------------------
HIV1$celltype_condition <- paste(
  HIV1$Marie_celltype,
  HIV1$condition,
  sep = " | "
)

Idents(HIV1) <- "celltype_condition"
DefaultAssay(HIV1) <- "RNA"

# Order groups: by cell type (same order as genes), then condition
ct_order     <- unique(top5_HIV$celltype)
cond_order   <- sort(unique(as.character(HIV1$condition)))
group_levels <- as.vector(t(outer(ct_order, cond_order, paste, sep = " | ")))
group_levels <- group_levels[group_levels %in% unique(as.character(Idents(HIV1)))]

# ---- 4. Compute average / percent expression (no DotPlot needed) ------------
# Seurat v5 uses `layer`, v4 uses `slot`
expr <- tryCatch(
  FetchData(HIV1, vars = genes_to_plot, layer = "data"),
  error = function(e) FetchData(HIV1, vars = genes_to_plot, slot = "data")
)

grp <- as.character(Idents(HIV1))
names(grp) <- colnames(HIV1)

df <- expr %>%
  rownames_to_column("cell") %>%
  mutate(group = grp[cell]) %>%
  filter(group %in% group_levels) %>%
  pivot_longer(-c(cell, group), names_to = "features.plot", values_to = "expr") %>%
  group_by(group, features.plot) %>%
  summarise(
    avg.exp = mean(expm1(expr)),
    pct.exp = 100 * mean(expr > 0),
    .groups = "drop"
  ) %>%
  group_by(features.plot) %>%
  mutate(avg.exp.scaled = as.numeric(scale(log1p(avg.exp)))) %>%
  ungroup() %>%
  rename(id = group)

df$features.plot <- factor(df$features.plot, levels = rev(genes_to_plot))
df$id            <- factor(df$id, levels = group_levels)

# ---- 5. Plot ----------------------------------------------------------------
p <- ggplot(df, aes(x = id, y = features.plot,
                    size = pct.exp, colour = avg.exp.scaled)) +
  geom_point() +
  scale_colour_distiller(palette = "RdBu", direction = -1,
                         limits = c(-2, 2.5), oob = scales::squish,
                         name = "Avg. exp.") +
  scale_size(range = c(0, 6), limits = c(0, 100), breaks = c(0, 25, 50, 75),
             name = "Pct. exp.") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title       = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3)
  )

print(p)

p <- ggplot(df, aes(x = id, y = features.plot,
                    size = pct.exp, colour = avg.exp.scaled)) +
  geom_point() +
  scale_colour_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1.5), oob = scales::squish,
    breaks = c(-1, 0, 1, 1.5), labels = c("≤ -1", "0", "1", "≥ 1.5"),
    name = "Avg. exp."
  ) +
  scale_size(range = c(0, 6), limits = c(0, 100), breaks = c(0, 25, 50, 75),
             name = "Pct. exp.") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3))


p
ggsave("E:/top2_HIV vs Mock_dotplot5.png", p,
       width  = max(8, 0.28 * length(group_levels) + 3),
       height = max(6, 0.22 * length(genes_to_plot) + 2))

###################################HIV vs Mock DOT PLOT #####################

library(Seurat); library(dplyr); library(ggplot2)

hm <- subset(HIV1, subset = condition %in% c("HIV", "Mock"))
Idents(hm) <- "Marie_celltype"

exclude <- "^MT-|^MTRNR|^RP[SL][0-9]|^RPLP|^SNHG|^UTY$|^USP9Y$|^PRKY$|^DDX3Y$|^KDM5D$|^EIF1AY$|^RPS4Y1$|^XIST$|^TSIX$"

# HIV vs Mock within each cell type (skip types with < 20 cells in either condition)
de <- lapply(levels(Idents(hm)), function(ct) {
  cells <- WhichCells(hm, idents = ct)
  if (sum(hm$condition[cells] == "HIV") < 20 || sum(hm$condition[cells] == "Mock") < 20) return(NULL)
  FindMarkers(hm, ident.1 = "HIV", ident.2 = "Mock", group.by = "condition",
              subset.ident = ct, min.pct = 0.1, logfc.threshold = 0.25) %>%
    tibble::rownames_to_column("gene") %>% mutate(celltype = ct)
}) %>% bind_rows() %>%
  filter(!grepl(exclude, gene, ignore.case = TRUE))

# How many genes change per cell type?
de %>% filter(p_val_adj < 0.05) %>% count(celltype, sign(avg_log2FC))

# Top 2 genes up in HIV per cell type
genes <- de %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  group_by(celltype) %>% slice_max(avg_log2FC, n = 2, with_ties = FALSE) %>%
  ungroup() %>% pull(gene) %>% unique()

# Keep only cell types that were tested (this removes the NA column)
tested <- unique(de$celltype)
hm2 <- subset(hm, subset = Marie_celltype %in% tested)
hm2$celltype_condition <- paste(hm2$Marie_celltype, hm2$condition, sep = " | ")

lv <- as.vector(t(outer(tested, c("Mock", "HIV"), paste, sep = " | ")))
lv <- lv[lv %in% unique(hm2$celltype_condition)]

# Dot plot data
df <- DotPlot(hm2, features = genes, group.by = "celltype_condition")$data
df$features.plot <- factor(df$features.plot, levels = rev(genes))
df$id <- factor(df$id, levels = lv)
df <- df %>% filter(!is.na(id))

# Plot
p <- ggplot(df, aes(id, features.plot, size = pct.exp, colour = avg.exp.scaled)) +
  geom_point() +
  scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1.5), oob = scales::squish,
                         breaks = c(-1, 0, 1, 1.5), labels = c("≤ -1", "0", "1", "≥ 1.5"),
                         name = "Avg. exp.") +
  scale_size(range = c(0, 6), limits = c(0, 100), breaks = c(0, 25, 50, 75),
             name = "Pct. exp.") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3))

print(p)
ggsave("E:/DE_HIV_vs_Mock_dotplot.png", p, width = 12, height = 8, dpi = 300)


avg <- AverageExpression(HIV1, features = genes_to_plot,
                         group.by = c("Marie_celltype", "condition"),
                         assays = "RNA", layer = "data")$RNA   # Seurat v5; use slot = "data" in v4

round(avg, 2)
range(df$avg.exp)             # real mean expression (non-log)
range(df$avg.exp.scaled)      # z-scores that set the dot colour
summary(df$avg.exp.scaled)

# Range per gene
df %>%
  group_by(features.plot) %>%
  summarise(min = min(avg.exp), max = max(avg.exp),
            min_z = min(avg.exp.scaled), max_z = max(avg.exp.scaled))
ggsave("E:/top2_HIV vs Mock_dotplot4.png", p,
       width  = max(8, 0.28 * length(group_levels) + 3),
       height = max(6, 0.22 * length(genes_to_plot) + 2))


genes_to_plot <- unique(top5_HIV$gene)

HIV1$celltype_condition <- paste(HIV1$Marie_celltype, HIV1$condition, sep = " | ")
Idents(HIV1) <- "celltype_condition"
DefaultAssay(HIV1) <- "RNA"

DotPlot(HIV1, features = genes_to_plot) +
  coord_flip() +
  scale_colour_distiller(palette = "RdBu", direction = -1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank())


library(Seurat); library(dplyr); library(ggplot2); library(RColorBrewer)

Idents(HIV1) <- "celltype_condition"          # your 19 cluster labels
DefaultAssay(HIV1) <- "RNA"

# 2 markers per cluster, genes ordered by cluster order
markers <- FindAllMarkers(HIV1, only.pos = TRUE)   # or your existing markers table

genes <- markers %>%
  filter(p_val_adj < 0.05) %>%
  mutate(cluster = factor(cluster, levels = levels(Idents(HIV1)))) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 2, with_ties = FALSE) %>%
  arrange(cluster) %>%
  pull(gene) %>% unique()

DotPlot(HIV1, features = genes, dot.scale = 6) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, "Spectral")),
                         name = "Exp.avg") +
  labs(x = "Features", y = "Identity", size = "Exp %") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
library(Seurat)
library(ggplot2)

# ============================================================
# Cluster annotations
# ============================================================

library(Seurat)
library(ggplot2)

# Define cluster annotations based on marker genes
cell_type <- c(
  "0"  = "Cytotoxic/Exhausted CD8+ T",
  "1"  = "CD4+ Helper T (Th2)",
  "2"  = "B cells",
  "3"  = "Naive/Memory CD4+ T",
  "4"  = "Effector Memory T",
  "5"  = "Gamma-Delta T / NK-T",
  "6"  = "Monocytes/Macrophages",
  "7"  = "Plasmacytoid DCs (pDC)",
  "8"  = "Cycling/Proliferating cells",
  "9"  = "High-Ribosomal/Metabolic cells",
  "10" = "NK cells",
  "11" = "Regulatory/Mature DCs",
  "12" = "Megakaryocytes/Platelets",
  "13" = "Stromal/Neuronal-like cells",
  "14" = "Granulocytes (Eosinophils/Basophils)",
  "15" = "Conventional DCs (cDC2)"
)


# ============================================================
# Make sure identities are the original Seurat clusters
# ============================================================

Idents(multiome) <- "seurat_clusters"
multiome<- RenameIdents(multiome, cell_type)

# Save annotations to object metadata for easy access
multiome$cell_type <- Idents(multiome)

# Visualize annotated clusters on UMAP
DimPlot(multiome, reduction = "wnn.umap", label = TRUE, pt.size = 0.5, label.size = 3.5) + 
  ggtitle("Cell-Type Annotations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# ============================================================
# Check cluster IDs
# ============================================================


# ============================================================
# Create cell-type metadata
# ============================================================

multiome$cell_type1 <- unname(
  cluster_annotations[as.character(multiome$seurat_clusters)]
)


# ============================================================
# Check annotations
# ============================================================

print("Cell type counts:")
print(table(multiome$cell_type1, useNA = "ifany"))


# ============================================================
# Set annotated cell types as active identities
# ============================================================

Idents(multiome) <- "cell_type"


# ============================================================
# Check final identities
# ============================================================

print("Final cell-type identities:")
print(levels(Idents(multiome)))

table(multiome$cell_type)
# ============================================================
# Visualize annotated WNN UMAP
# ============================================================

DimPlot(
  multiome,
  reduction = "wnn.umap",
  group.by = "cell_type", split.by="condition",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Annotated Cell Types") +
  theme_classic()


# ============================================================
# Optional: save the annotated Seurat object
# ============================================================

 saveRDS(
   multiome,
   file = "F:/Elements/All_3_Donors_multiome_annotated.rds"
 )



library(Seurat)
library(ggplot2)

library(Seurat)
library(ggplot2)

# ============================================================
# 1. Cell-cycle genes
# ============================================================

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes


# ============================================================
# 2. RNA Assay
# ============================================================

DefaultAssay(multiome) <- "RNA"


# ============================================================
# 3. Cell-cycle scoring
# ============================================================

multiome <- CellCycleScoring(
  object = multiome,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE
)


# ============================================================
# 4. Check cell-cycle information
# ============================================================

head(
  multiome@meta.data[
    ,
    c("S.Score", "G2M.Score", "Phase")
  ]
)

table(multiome$Phase)


# ============================================================
# 5. Visualize cell-cycle phase on WNN UMAP
# ============================================================

DimPlot(
  multiome,
  reduction = "wnn.umap",
  group.by = "Phase", split.by="condition",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Cell Cycle Phase") +
  theme_classic()


# ============================================================
# 6. Visualize S-phase score
# ============================================================

FeaturePlot(
  multiome,
  features = "S.Score",
  reduction = "wnn.umap"
) +
  ggtitle("S Phase Score") +
  theme_classic()


# ============================================================
# 7. Visualize G2/M score
# ============================================================

FeaturePlot(
  multiome,
  features = "G2M.Score",
  reduction = "wnn.umap"
) +
  ggtitle("G2/M Phase Score") +
  theme_classic()


#
DimPlot(
  multiome,
  reduction = "wnn.umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Annotated Cell Types") +
  theme_classic()


# ============================================================
# 14. Cell cycle phase by cell type
# ============================================================

DimPlot(
  multiome,
  reduction = "wnn.umap",
  group.by = "Phase",
  split.by = "cell_type"
)


DimPlot(
  multiome,
  reduction = "wnn.umap",
  group.by = "Phase"
)

# ============================================================
# 15. Cell-cycle phase distribution by cell type
# ============================================================

cellcycle_table <- table(
  multiome$cell_type,
  multiome$Phase
)

print(cellcycle_table)


# ============================================================
# 16. Percentage of each cell type in each cell-cycle phase
# ============================================================

cellcycle_percent <- prop.table(
  cellcycle_table,
  margin = 1
) * 100

print(round(cellcycle_percent, 2))


# ============================================================
# 17. Violin plots for S and G2/M scores
# ============================================================

VlnPlot(
  multiome,
  features = c("S.Score", "G2M.Score"),
  group.by = "cell_type",
  pt.size = 0
) 


Assays(multiome)
DefaultAssay(multiome) <- "ATAC"
DefaultAssay(multiome) <- "ATAC"

library(Signac)
library(Seurat)

# Plot ATAC accessibility across cell cycle phases at the MKI67 locus
CoveragePlot(
  object = multiome,
  region = "MX1",
  group.by = "condition",
  
  # Groups cells by G1, S, G2M cell cycle phase
  annotation = TRUE,           # Displays gene structure
  peaks = TRUE,                # Displays called peak regions
  extend.upstream = 2000,
  extend.downstream = 2000
)

CoveragePlot(
  object = multiome,
  region = "MX1",
  group.by ="seurat_clusters",
  
  annotation = TRUE,           # Displays gene structure
  peaks = TRUE,                # Displays called peak regions
  extend.upstream = 2000,
  extend.downstream = 2000
)
CoveragePlot(
  multiome,
  region = "TXN",
  group.by = "condition",
  extend.upstream = 3000,
  extend.downstream = 3000
)



library(Seurat)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
install.packages("EnhancedVolcano")
BiocManager::install("EnhancedVolcano")
# Ensure active assay is RNA or SCT
DefaultAssay(multiome) <- "RNA"

# -----------------------------------------------------------------------------
# 1. REMOVE RIBOSOMAL GENES
# -----------------------------------------------------------------------------
# Identify ribosomal genes (RPL and RPS families)
ribo_genes <- grep(pattern = "^RP[SL]", x = rownames(multiome), value = TRUE)
cat("Found", length(ribo_genes), "ribosomal genes to remove.\n")

# Filter out ribosomal genes from the feature matrix
non_ribo_genes <- setdiff(rownames(multiome), ribo_genes)
multiome_filtered <- subset(multiome, features = non_ribo_genes)

# -----------------------------------------------------------------------------
# 2. SETUP OUTPUT DIRECTORIES
# -----------------------------------------------------------------------------
dir.create("F:/Elements/DEG_Results", showWarnings = FALSE)
dir.create("F:/Elements/Volcano_Plots", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 3. CLUSTER-WISE DEG ANALYSIS (HIV vs. Mock) & VOLCANO PLOTS
# -----------------------------------------------------------------------------
# Set active identity to cluster column (e.g., 'seurat_clusters' or 'major_lineage')
library(Seurat)
library(EnhancedVolcano)
library(ggplot2)

# ---------------------------------------------------------
# SETTINGS
# ---------------------------------------------------------

# Make sure cluster identities are set
Idents(multiome_filtered) <- "seurat_clusters"



# Process clusters 0-18
clusters <- as.character(0:15)


# ---------------------------------------------------------
# LOOP THROUGH ALL CLUSTERS
# ---------------------------------------------------------

for (cls in clusters) {
  
  cat("\n========================================\n")
  cat("Processing Cluster:", cls, "\n")
  cat("========================================\n")
  
  # -------------------------------------------------------
  # SUBSET CLUSTER
  # -------------------------------------------------------
  
  sub_obj <- subset(
    multiome_filtered,
    idents = cls
  )
  
  cat("Total cells:", ncol(sub_obj), "\n")
  
  
  # -------------------------------------------------------
  # CHECK HIV / MOCK
  # -------------------------------------------------------
  
  cond_counts <- table(sub_obj$condition)
  
  cat("Condition counts:\n")
  print(cond_counts)
  
  # Check that both conditions exist
  if (!all(c("HIV", "ssDNA") %in% names(cond_counts))) {
    
    cat(
      "Skipping cluster", cls,
      "- HIV or Mock is missing.\n"
    )
    
    next
  }
  
  # Require at least 3 cells per condition
  if (min(cond_counts[c("HIV", "ssDNA")]) < 3) {
    
    cat(
      "Skipping cluster", cls,
      "- fewer than 3 cells in HIV or Mock.\n"
    )
    
    next
  }
  
  
  # -------------------------------------------------------
  # CHECK DONOR
  # -------------------------------------------------------
  
  if (!"donor" %in% colnames(sub_obj@meta.data)) {
    
    cat(
      "Skipping cluster", cls,
      "- donor column not found.\n"
    )
    
    next
  }
  
  
  # -------------------------------------------------------
  # SET CONDITION AS IDENTITY
  # -------------------------------------------------------
  
  Idents(sub_obj) <- "condition"
  
  
  # -------------------------------------------------------
  # HIV vs MOCK DEG
  # -------------------------------------------------------
  
  cat("Running HIV vs Mock DE analysis...\n")
  
  deg_results <- FindMarkers(
    object = sub_obj,
    ident.1 = "HIV",
    ident.2 = "ssDNA",
    test.use = "LR",
    latent.vars = "donor",
    min.pct = 0.1,
    logfc.threshold = 0.25
  )
  
  
  # -------------------------------------------------------
  # ADD GENE NAMES
  # -------------------------------------------------------
  
  deg_results$gene <- rownames(deg_results)
  
  
  # -------------------------------------------------------
  # REMOVE RIBOSOMAL GENES
  # -------------------------------------------------------
  
  ribosomal_genes <- grepl(
    "^RPL|^RPS|^MRPL|^MRPS",
    deg_results$gene,
    ignore.case = TRUE
  )
  
  deg_results <- deg_results[!ribosomal_genes, ]
  
  
  # -------------------------------------------------------
  # SORT BY ADJUSTED P VALUE
  # -------------------------------------------------------
  
  deg_results <- deg_results[
    order(deg_results$p_val_adj),
  ]
  
  
  # -------------------------------------------------------
  # SAVE DEG RESULTS
  # -------------------------------------------------------
  
  csv_filename <- paste0(
    "F:/Elements/DEG_Results/",
    "DEG_Cluster_",
    cls,
    "_HIV_vs_ssDNA.csv"
  )
  
  write.csv(
    deg_results,
    file = csv_filename,
    row.names = FALSE
  )
  
  cat(
    "Saved DEG table:",
    csv_filename,
    "\n"
  )
  
  
  # -------------------------------------------------------
  # PREPARE VOLCANO DATA
  # -------------------------------------------------------
  
  volcano_data <- deg_results[
    !is.na(deg_results$p_val_adj) &
      !is.na(deg_results$avg_log2FC) &
      deg_results$p_val_adj > 0,
  ]
  
  
  # -------------------------------------------------------
  # CHECK WHETHER THERE ARE GENES TO PLOT
  # -------------------------------------------------------
  
  if (nrow(volcano_data) == 0) {
    
    cat(
      "No valid genes for volcano plot in cluster",
      cls,
      "\n"
    )
    
    next
  }
  
  
  # -------------------------------------------------------
  # VOLCANO PLOT
  # -------------------------------------------------------
  
  volcano_p <- EnhancedVolcano(
    toptable = volcano_data,
    
    lab = volcano_data$gene,
    
    x = "avg_log2FC",
    
    y = "p_val_adj",
    
    pCutoff = 0.05,
    
    FCcutoff = 0.5,
    
    pointSize = 2.0,
    
    labSize = 4.0,
    
    title = paste(
      "Cluster",
      cls,
      ": HIV vs ssDNA"
    ),
    
    subtitle = paste(
      "Ribosomal genes excluded |",
      "Controlled for donor"
    ),
    
    caption = paste(
      "log2FC cutoff = 0.5 |",
      "FDR cutoff = 0.05"
    ),
    
    legendPosition = "right"
  )
  
  
  # -------------------------------------------------------
  # SAVE VOLCANO PLOT
  # -------------------------------------------------------
  
  png_filename <- paste0(
    "F:/Elements/Volcano_Plots/",
    "Volcano_Cluster_",
    cls,
    "_HIV_vs_ssDNA.png"
  )
  
  ggsave(
    filename = png_filename,
    plot = volcano_p,
    width = 8,
    height = 7,
    dpi = 300
  )
  
  cat(
    "Saved volcano plot:",
    png_filename,
    "\n"
  )
}


# ---------------------------------------------------------
# FINISHED
# ---------------------------------------------------------

cat("\n========================================\n")
cat("ALL CLUSTERS PROCESSED\n")
cat("========================================\n")
}

cat("\n--- All Cluster DEG Analysis & Volcano Plots Completed! ---\n")


Assays(multiome)
Idents(multiome) <- "seurat_clusters"

# Subset cluster 1
cluster0 <- subset(
  multiome,
  idents = "0"
)

cluster1
colnames(cluster0@meta.data)
table(cluster0$condition)
cluster0_HIV_Mock <- subset(
  cluster0,
  subset = condition %in% c("HIV", "Mock")
)
Assays(cluster1_HIV_Mock)
table(cluster1_HIV_Mock$condition)
cluster0_HIV_Mock 
CoveragePlot(
  object = cluster0_HIV_Mock,
  region = "IFI44L",
  assay = "ATAC",
  group.by = "condition",
  extend.upstream = 3000,
  extend.downstream = 3000
)

########################Linkplot###################################3  
  
  DefaultAssay(multiome) <- "ATAC"
  
  multiome <- RegionStats(
    object = multiome,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )
  multiome <- LinkPeaks(
    object = multiome,
    peak.assay = "ATAC",
    expression.assay = "RNA"
  )
  
  Links(multiome)
  
  CoveragePlot(
    object = cluster0_HIV_Mock,
    assay = "ATAC",
    region = "MX1",
    group.by = "condition",
    extend.upstream = 3000,
    extend.downstream = 3000,
    links = TRUE
  )
  saveRDS(
    multiome,
    file = "E:/All_4_Donors_multiome_linkpeaks.rds"
  )
  cluster9 <- subset(
    multiome,
    idents = "9"
  )
  

  cluster9_HIV_Mock <- subset(
    cluster9,
    subset = condition %in% c("HIV", "Mock"))
  
  
  CoveragePlot(
    object = cluster9_HIV_Mock,
    assay = "ATAC",
    region = "CD52",
    group.by = "condition",
    extend.upstream = 3000,
    extend.downstream = 3000,
    links = TRUE
  )
  # Calculate peak co-accessibility (e.g., using Cicero or Signac's LinkPeaks)
  cluster9_HIV_Mock <- LinkPeaks(
    object = cluster9_HIV_Mock,
    peak.assay = "ATAC",
    expression.assay = "RNA", # Optional: link peaks to gene expression if RNA data exists
    genes.use = "GZMB",
    distance = 500000
  )

linkplot <- CoveragePlot(
    object = cluster9_HIV_Mock,
    assay = "ATAC",
    region = "GZMB",
    group.by = "condition",
    extend.upstream = 30000,
    extend.downstream = 30000,
    links = TRUE
  )
  
  # Display the plot
  print(linkplot)
  
  
  
  cluster12 <- subset(
    multiome,
    idents = "12"
  )
  
  
  cluster12_HIV_Mock <- subset(
    cluster12,
    subset = condition %in% c("HIV", "Mock"))
  
  
  CoveragePlot(
    object = cluster12_HIV_Mock,
    assay = "ATAC",
    region = "TMSB10",
    group.by = "condition",
    extend.upstream = 3000,
    extend.downstream = 3000,
    links = TRUE
  )
  # Calculate peak co-accessibility (e.g., using Cicero or Signac's LinkPeaks)
  cluster10_HIV_Mock <- LinkPeaks(
    object = cluster10_HIV_Mock,
    peak.assay = "ATAC",
    expression.assay = "RNA", # Optional: link peaks to gene expression if RNA data exists
    genes.use = "IFI27",
    distance = 5000
  )
  
  linkplot <- CoveragePlot(
    object = cluster9_HIV_Mock,
    assay = "ATAC",
    region = "GZMB",
    group.by = "condition",
    extend.upstream = 30000,
    extend.downstream = 30000,
    links = TRUE
  )
  
  # Display the plot
  print(linkplot)
  
  
  library(JASPAR2020)
  library(TFBSTools)
  library(motifmatchr)
  
  # 1. Get vertebrate JASPAR transcription factor motifs
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", species = 9606, all_versions = FALSE)
  )
  
  # 2. Add motif annotations to your Seurat/Signac object
  cluster0_HIV_Mock <- AddMotifs(
    object = cluster0_HIV_Mock,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    pfm = pfm
  )

##############################################DA Peaks #########################################################################  
  library(Seurat)
  library(Signac)
  library(dplyr)
  
  HIV1<-readRDS("E:/Allclusters_Marie_celltype.rds")
  HIV1
  DimPlot(
    HIV1,
    reduction = "wnn.umap",
    group.by = "Marie_celltype",
    label = TRUE
  )
  DimPlot(
    HIV1,
    reduction = "wnn.umap",
    group.by = "Marie_celltype",
    #label = TRUE,
    repel = TRUE,
    label.size = 5
  ) +
    ggtitle("WNN UMAP by Cell Type") +
    theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14)
    )
  colnames(HIV@meta.data)
  library(Seurat)
  library(ggplot2)
  
  p1 <- DimPlot(
    HIV1,
    reduction = "wnn.umap",
    group.by = "Marie_celltype",
    #label = TRUE,
    repel = TRUE,
    label.size = 4,
    pt.size = 0.3
  ) +
    labs(title = "Integrated RNA–ATAC Cell Atlas",
         x = "WNN UMAP 1",
         y = "WNN UMAP 2",
         color = "Cell type") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(
        size = 18, face = "bold", hjust = 0.5
      ),
      axis.text = element_text(size = 11),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10)
    )
  
  p1
  
  ggsave("E:/WNN_celltype_UMAP.pdf", p1,
         width = 10, height = 7)
  ggsave("WNN_celltype_UMAP.png", p1,
         width = 10, height = 7, dpi = 600)
  
 ########################Cell Counts Plot######################## 
  library(dplyr)
  library(ggplot2)
  library(scales)
  
  df <- HIV1@meta.data %>%
    count(condition, seurat_clusters, name = "cells") %>%
    group_by(condition) %>%
    mutate(proportion = cells / sum(cells)) %>%
    ungroup()
  
  p2 <- ggplot(
    df,
    aes(x = condition, y = proportion,
        fill = seurat_clusters)
  ) +
    geom_col(width = 0.75, color = "white",
             linewidth = 0.15) +
    scale_y_continuous(
      labels = percent_format(),
      expand = c(0, 0)
    ) +
    labs(
      title = "Cell-Type Composition by Condition",
      x = NULL,
      y = "Proportion of cells",
      fill = "Cell type"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(
        size = 18, face = "bold", hjust = 0.5
      ),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  p2
  paper_colors <- c(
    "#0072B2", # Blue
    "#E69F00", # Orange
    "#009E73", # Green
    "#D55E00", # Vermillion
    "#CC79A7", # Pink-purple
    "#56B4E9", # Sky blue
    "#F0E442", # Yellow
    "#000000", # Black
    "#999999", # Grey
    "#332288", # Indigo
    "#88CCEE", # Light blue
    "#44AA99", # Teal
    "#AA4499", # Purple
    "#DDCC77", # Sand
    "#882255", # Wine
    "#117733", # Dark green
    "#661100", # Brown
    "#6699CC" , # Steel blue
    "#787655",
    "#993344"
  )
  
  p2 <- p2 +
    scale_fill_manual(
      values = setNames(
        paper_colors,
        levels(factor(HIV1$seurat_clusters))
      )
    )
  
  p2
  
  library(ggplot2)
  library(dplyr)
  
  cell_counts <- HIV1@meta.data %>%
    filter(!is.na(Marie_celltype)) %>%
    count(Marie_celltype, name = "Cell_count") %>%
    arrange(Cell_count)
  
  p <- ggplot(
    cell_counts,
    aes(x = Cell_count, y = reorder(Marie_celltype, Cell_count))
  ) +
    geom_col(fill = "#2878B5", width = 0.75) +
    geom_text(
      aes(label = Cell_count),
      hjust = -0.15,
      size = 3.5
    ) +
    labs(
      title = "Cell Number by Cell Type",
      x = "Number of cells",
      y = "Cell type"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(
        size = 17, face = "bold", hjust = 0.5
      ),
      axis.text.y = element_text(size = 10)
    )
  
  p
  
  ggsave("E:/Cell_number_by_celltype.pdf", p,
         width = 10, height = 7)
  
  ggsave("E:/Celltype_composition.pdf", p2,
         width = 9, height = 6)
  
  
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  
p_cluster <- ggplot(
  cluster_info,
  aes(
    x = Cells,
    y = reorder(Cluster, Cells),
    fill = Marie_celltype
  )
) +
  geom_col(width = 0.75, show.legend = FALSE) +
  geom_text(
    aes(label = Cells),
    hjust = -0.15,
    size = 3.5
  ) +
  labs(
    title = "WNN Cluster Numbers and Cell Types",
    x = "Number of cells",
    y = "Cluster — Cell type"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      size = 17, face = "bold", hjust = 0.5
    ),
    axis.text.y = element_text(size = 10)
  )

p_cluster

ggsave(
  "E:/WNN_cluster_celltype.pdf",
  p_cluster,
  width = 12,
  height = 8
)
p_cluster <- ggplot(
  cluster_info,
  aes(
    x = Cells,
    y = reorder(Cluster, Cells),
    fill = Marie_celltype
  )
) +
  geom_col(width = 0.75, show.legend = FALSE) +
  geom_text(
    aes(label = Cells),
    hjust = -0.15,
    size = 3.5
  ) +
  labs(
    title = "WNN Cluster Numbers and Cell Types",
    x = "Number of cells",
    y = "Cluster — Cell type"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      size = 17, face = "bold", hjust = 0.5
    ),
    axis.text.y = element_text(size = 10)
  )

p_cluster

ggsave(
  "E:/WNN_cluster_celltype.pdf",
  p_cluster,
  width = 12,
  height = 8
)
  
  ggsave("E:/WNN_cluster_celltype.pdf",
         p_cluster, width = 12, height = 7)
  
  table(HIV$condition)
  # ============================================================
  # ATAC differential accessibility: HIV vs Mock in cluster 0
  # ============================================================
  
  # 1. Check the condition labels
  table(HIV$condition,HIV$seurat_clusters)
  
  # 2. Check cluster labels
  table(HIV$seurat_clusters)
  HIV<-HIV1
  # ------------------------------------------------------------
  # 3. Subset cluster 0
  # ------------------------------------------------------------
  Idents(HIV)<-"seurat_clusters"
  HIV_cluster0 <- subset(
    HIV,
    subset = seurat_clusters == 4
  )
  
  # Check number of cells
  table(HIV_cluster0$condition)
  
  # ------------------------------------------------------------
  # 4. Set ATAC as the active assay
  # ------------------------------------------------------------
  
  DefaultAssay(HIV_cluster0) <- "ATAC"
  
  # Set identities to HIV / Mock
  Idents(HIV_cluster0) <- HIV_cluster0$condition
  
  # Check identities
  table(Idents(HIV_cluster0))
  
  # ------------------------------------------------------------
  # 5. Differentially accessible peaks
  #    HIV vs Mock
  # ------------------------------------------------------------
  
  DA_HIV_vs_Mock_cluster0 <- FindMarkers(
    HIV_cluster0,
    ident.1 = "HIV",
    ident.2 = "Mock",
    only.pos = TRUE,
    test.use = "LR"
  )
  
  # ------------------------------------------------------------
  # 6. Sort by adjusted p-value
  # ------------------------------------------------------------
  
  DA_HIV_vs_Mock_cluster0 <- DA_HIV_vs_Mock_cluster0[
    order(DA_HIV_vs_Mock_cluster0$p_val_adj),
  ]
  
  # View top 20 HIV-enriched peaks
  head(DA_HIV_vs_Mock_cluster0, 20)
  
  # ------------------------------------------------------------
  # 7. Significant peaks
  # ------------------------------------------------------------
  
  sig_peaks_cluster0 <- subset(
    DA_HIV_vs_Mock_cluster0,
    p_val_adj < 0.05 & avg_log2FC > 0.25
  )
  
  # Number of significant peaks
  nrow(sig_peaks_cluster0)
  
  # Top significant peaks
  head(sig_peaks_cluster0, 20)
  # ============================================================
  # 4. SELECT SIGNIFICANT HIV-ENRICHED PEAKS
  # ============================================================
  
  sig_peaks <- rownames(
    subset(
      DA_HIV_vs_Mock_cluster0,
      p_val_adj < 0.05 &
        avg_log2FC > 0.25
    )
  )
  
  length(sig_peaks)
  HIV_Mock_obj <- subset(HIV_cluster0, condition %in% c("HIV", "Mock"))
  # Top 20
  top_peaks <- sig_peaks[
    1:min(20, length(sig_peaks))
  ]
  
  top_peaks
  # ------------------------------------------------------------
  # 8. Save results
  # ------------------------------------------------------------
  
  write.csv(
    DA_HIV_vs_Mock_cluster0,
    "E:/Cluster3_ATAC_DA_HIV_vs_Mock.csv"
  )
  
  write.csv(
    sig_peaks_cluster0,
    "E:/Cluster3_ATAC_DA_HIV_enriched.csv"
  )
  library(Signac)
  library(Seurat)
  
  # 1. Annotate significant peaks to find the closest genes

  annotated_peaks <- ClosestFeature(HIV_cluster0, regions = rownames(sig_peaks_cluster0))
  
  # 2. Combine original differential statistics with the gene annotations
  sig_peaks_with_genes <- cbind(sig_peaks_cluster0, annotated_peaks)
  
  # 3. View the first few rows to check that gene names/symbols are now included
  head(sig_peaks_with_genes)
  
  # 4. Save this rich table to your E:/ drive
  write.csv(
    sig_peaks_with_genes,
    file = "E:/Cluster3_ATAC_DA_HIV_enriched_with_Genes.csv",
    row.names = TRUE
  )
  
  # ============================================================
  # 6. COVERAGE PLOT
  # ============================================================
  
  # Plot the strongest DA peak
  CoveragePlot(
    HIV_cluster0,
    region = top_peaks[1],
    assay = "ATAC",
    group.by = "condition"
  )
  
  
  # ============================================================
  # 7. COVERAGE PLOT FOR TOP 5 PEAKS
  # ============================================================
  
 CoveragePlot(
    HIV_cluster0,
    region = top_peaks[1:min(5, length(top_peaks))],
    assay = "ATAC",
    group.by = "condition"
  )
  p1 <- CoveragePlot(
    object = HIV_Mock_obj,
    region = top_peaks[1:min(5, length(top_peaks))],
    assay = "ATAC",
    group.by = "condition"
  )
  p1
  # 3. Save the clean plot
  ggsave(
    filename = "E:/p3.png",
    plot = p1,
    width = 10,
    height = 6,
    units = "in"
  )
  ggsave(
    filename = "E:/p1.png",
    plot = p1,
    width = 10,
    height = 6,
    units = "in"
  )
  
  # From scRNAseq
  genes_to_plot <- c("IFI44L", "IFI6", "TMSB10", "IFI44")
  genes_to_plot <- c("ISG15", "IFITM3")
  # Adjust based on your top hits
  genes_to_plot <- c("PLCG2", "IFI44L", "CCL3","TMSB10")
  genes_to_plot <- c("PRKCA","SKAP1")
  
  # Generate the multi-panel coverage plot
  p2<-CoveragePlot(
    object = HIV_Mock_obj,
    region = genes_to_plot,
    assay = "ATAC",
    group.by = "condition"
  )
  
  p2
  ggsave(
    filename = "E:/p2.png",
    plot = p2,
    width = 10,
    height = 6,
    units = "in"
  )
  # ============================================================
  # 8. CHECK PEAK-GENE LINKS
  # ============================================================
  
  Links(HIV_cluster0)
  
  
  # ============================================================
  # 9. FIND GENES ASSOCIATED WITH YOUR PEAKS
  # ============================================================
  
  # Check the gene annotations available in the ATAC assay
  Annotation(HIV_cluster0[["ATAC"]])
  
  
  # ============================================================
  # 10. LINK PLOT
  # ============================================================
  
  library(BSgenome.Hsapiens.UCSC.hg38)
  
  # 1. Compute region statistics for the peaks
  HIV_cluster0 <- RegionStats(HIV_cluster0, genome = BSgenome.Hsapiens.UCSC.hg38)
  
  # 2. #peak-to-gene regulatory links
  HIV_cluster0 <- LinkPeaks(
    object = HIV_cluster0,
    peak.assay = "ATAC",
    expression.assay = "RNA" # Make sure this matches your RNA assay name
  )
  
  ######################################################################
  # Saving
  #####################################################################
  # 1. Extract the links from the Seurat object as a data frame
  all_links_df <- as.data.frame(Links(HIV_cluster0))
  
  # 2. View the first few rows to check the columns (e.g., peak, gene, score, correlation)
  head(all_links_df)
  
  # 3. Save 
  write.csv(
    all_links_df,
    file = "E:/Cluster_3_ATAC_Peak_to_Gene_Links.csv",
    row.names = FALSE
  )
  
  
  # Filter the links to only include your specific key genes
  filtered_links <- all_links_df[all_links_df$gene %in% c("IFI44L", "IFI6", "ISG15", "IFITM3", "VIM"), ]
  
  # Save just the filtered list
  write.csv(filtered_links, file = "E:/Cluster_3_ISG_Peak_to_Gene_Links.csv", row.names = FALSE)
  # Create individual tracks
  # Check if your specific genes have any calculated links
  target_links <- all_links_df[all_links_df$gene %in% c("IFI44L", "IFI6", "ISG15", "IFITM3", "VIM"), ]
  print(target_links)
  # Subset to HIV vs Mock if needed
  HIV_Mock_obj <- subset(HIV_cluster0, condition %in% c("HIV", "Mock"))
  library(Signac)
  library(Seurat)
  print(head(Links(target_links)))
  # Define the exact genomic region covering IFI44L and its linked peak
  target_region <- "chr11-279039-279797"
  
  # Generate Coverage Plot with explicit coordinates
  cov_plot <- CoveragePlot(
    object = HIV_Mock_obj,
    region = target_region,
    assay = "ATAC",
    group.by = "condition"
  )
  
  # Generate Link Plot with the same explicit coordinates
  link_plot <- LinkPlot(
    object = HIV_cluster0,
    region = target_region,
    assay = "ATAC"
  )
  link_plot
  # Combine them
  final_plot <- CombineTracks(
    plotlist = list(cov_plot, link_plot),
    heights = c(4, 1)
  )
  final_plot
  # 6. Save using ggsave
  ggsave(
    filename = "E:/VIM_Coverage_and_Links.png",
    plot = final_plot,
    width = 10,
    height = 6,
    units = "in"
  )
  
  message("Plot successfully saved to E:/IFI44L_Coverage_and_Links.pdf")

  
  message("Link plots successfully saved to E:/Cluster_0_ISG_LinkPlots.pdf")
  cov_track <- CoveragePlot(object = HIV_Mock_obj, region = "IFI44L")
  link_track <- LinkPlot(object = HIV_Mock_obj, region = "IFI44L")
  
  # Combine them vertically
  CombineTracks(plotlist = list(cov_track, link_track))

  table(multiome$condition, multiome$donor)
  
####################################################################################################
# Micke's Suggestion
####################################################################################################
  Assays(HIV1)
  head(HIV1@meta.data)
  colnames(HIV@meta.data)
  Reductions(HIV)
  
  HIV[["lsi"]]
  
  DimPlot(
    HIV,
    reduction = "lsi",
    dims = c(2, 3),
    group.by = "cell_type"
  )
  
  library(Seurat)
  library(Signac)
  
  DefaultAssay(HIV1) <- "ATAC"
  
  # Skip LSI component 1 (usually correlates with sequencing depth).
  # Check with DepthCor(HIV, reduction = "lsi") and use harmony_atac if that's what you integrated on.
  HIV1 <- FindNeighbors(HIV1, reduction = "harmony_atac", dims = 2:30,
                       graph.name = c("atac_nn", "atac_snn"))
  HIV1 <- FindClusters(HIV1, graph.name = "atac_snn", algorithm = 3,  # SLM, standard for ATAC
                      resolution = 0.5, verbose = FALSE)
  HIV1$atac_clusters <- HIV1$seurat_clusters
  
  HIV1 <- RunUMAP(HIV1, reduction = "harmony_atac", dims = 2:30,
                 reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  
  DimPlot(HIV1, reduction = "umap.atac", group.by = "atac_clusters", label = TRUE)
  DimPlot(HIV1, reduction = "umap.atac", group.by = "condition")   # check conditions mix within clusters 
  
  # Option A: RNA markers per ATAC cluster
  DefaultAssay(HIV) <- "RNA"
  HIV <- NormalizeData(HIV)   # skip if already done
  Idents(HIV) <- "atac_clusters"
  markers <- FindAllMarkers(HIV, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.25)
  dim(markers)
  head(markers)
  table(markers$cluster)   # all 20 clusters should appear
  
  # RNA clusters: use whichever column holds them (check colnames(HIV@meta.data))
  table(HIV1$seurat_clusters)
  table(HIV$atac_clusters)
  tab <- table(ATAC = HIV$atac_clusters, RNA = HIV1$seurat_clusters)
  round(prop.table(tab, 1), 2)   # each ATAC cluster (row): which RNA clusters it falls into
 
  # Heatmap if you want to see it
  pheatmap::pheatmap(prop.table(tab, 1), cluster_rows = TRUE, cluster_cols = TRUE)
  
  
  DotPlot(HIV, features = c("CD3E","CD4","CD8A","IL7R","CCR7","NKG7","GNLY",
                            "CD14","LYZ","FCGR3A","MS4A1","CD79A","PPBP"),
          group.by = "atac_clusters") + RotatedAxis()
  
  library(dplyr)
  
  top <- markers %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 10) %>%
    summarise(top_genes = paste(gene, collapse = ", "))
  
  print(top, n = 20, width = Inf)
  
  DefaultAssay(HIV) <- "RNA"
  DotPlot(HIV1, features = c("CD3E","CD4","CD8A","IL7R","CCR7","SELL","GZMK","GZMH",
                            "NKG7","GNLY","NCAM1","TRDC","CD14","LYZ","FCGR3A",
                            "MS4A1","CD79A","PPBP","HBB","MKI67"),
          group.by = "atac_clusters") + RotatedAxis()
   DotPlot(HIV1, group.by = "atac_clusters",
               features = c(
                 # cDC1
                 "CLEC9A", "XCR1", "CADM1",
                 # cDC2
                 "CD1C", "CLEC10A", "FCER1A",
                 # Activated / Migratory DCs
                 "LAMP3", "CCR7", "CCL19", "IDO1",
                 # Plasmacytoid DCs (pDC)
                 "IL3RA", "LILRA4", "TCF4", "IRF8",
                 # Classical / Non-classical Monocytes
                 "CD14", "S100A8", "FCGR3A", "HLA-DRA"
               )) + 
    RotatedAxis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
  
  ggsave("results/dotplot_comprehensive_DC.pdf", p, width = 12, height = 6)
  # QC for the odd clusters (10, 14, 8, 19, 18)
  VlnPlot(HIV, c("nCount_RNA","nFeature_RNA","nCount_ATAC","percent.mt"),
          group.by = "atac_clusters", pt.size = 0)
  
  # Median QC per cluster
  library(dplyr)
  HIV@meta.data %>%
    group_by(atac_clusters) %>%
    summarise(n = n(),
              rna_umi = median(nCount_RNA),
              genes   = median(nFeature_RNA),
              atac    = median(nCount_ATAC),
              mt      = median(percent.mt)) %>%
    print(n = 25)
  
  # Doublet scoring (RNA-based, run per sample/library if you have that column)
  library(scDblFinder)
  sce <- as.SingleCellExperiment(DietSeurat(HIV, assays = "RNA"))
  sce <- scDblFinder(sce, samples = HIV$sample)   # use your library/sample column
  HIV$dbl_score <- sce$scDblFinder.score
  HIV$dbl_class <- sce$scDblFinder.class
  
  table(HIV$atac_clusters, HIV$dbl_class)
  VlnPlot(HIV, "dbl_score", group.by = "atac_clusters", pt.size = 0)
  ###############################################################################
  # Cluster on ATAC only -> label -> within-cluster RNA DE between conditions
  # -> test whether DE genes have "primed" (accessible, non-differential) chromatin
  #
  # Starting point: Seurat/Signac object `HIV` with assays RNA + ATAC,
  # `atac_clusters` already computed, `dbl_score` / `dbl_class` from scDblFinder.
  ###############################################################################
  
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(Matrix)
  library(DESeq2)
  library(GenomicRanges)
  library(ggplot2)
  
  ###############################################################################
  # Cluster on ATAC only -> label -> within-cluster RNA DE between conditions
  # -> test whether induced genes have "primed" (open, unchanged) promoters
  #
  # v2: fixes for sparse.model.matrix, masked first(), raw vs shrunken LFC,
  #     optional ashr, tibble printing, stricter primed definition, diagnostics.
  #
  # Starting point: Seurat/Signac object `HIV` with RNA + ATAC assays,
  # `atac_clusters`, `dbl_class` (scDblFinder) already in the metadata.
  ###############################################################################
  
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(DESeq2)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)          # load LAST so dplyr verbs win over S4Vectors/GenomicRanges
  
  table(HIV1$condition)
  # ----------------------------- 0. SETTINGS (EDIT) ----------------------------
  SAMPLE_COL <- "sample"        # one value per biological sample/library
  DONOR_COL  <- NULL            # e.g. "donor" if paired design, else NULL
  COND_COL   <- "condition"
  COND_CASE  <- "HIV"      # <- must match a real level of HIV$condition
  COND_CTRL  <- "Mock"          # <- must match a real level of HIV$condition
  MIN_CELLS  <- 20              # min cells per (celltype x sample) pseudobulk
  #DROP_CLUSTERS <- c("8", "9", "10", "18", "20")
  
  # Primed-chromatin thresholds
  PADJ_RNA     <- 0.1           # exploratory DESeq2 cutoff
  LFC_RNA      <- 0.5           # raw (unshrunk) log2FC
  MIN_BASEMEAN <- 10
  PROM_UP      <- 2000          # promoter window upstream of TSS
  PROM_DOWN    <- 500
  PROTEIN_CODING_ONLY <- TRUE
  OPEN_PCT     <- 0.05          # promoter peak open in >= 5% of cells in each condition
  STABLE_LFC   <- 0.25          # |ATAC log2FC| below this = "stable"
  CPM_SILENT   <- 10            # control CPM below this = "silent at baseline"
  
  dir.create("E:/results", showWarnings = FALSE)
  
  stopifnot(all(c(SAMPLE_COL, COND_COL, "atac_clusters") %in%
                  colnames(HIV1@meta.data)))
  cat("Condition levels found:", paste(unique(HIV1@meta.data[[COND_COL]]), collapse = ", "), "\n")
  stopifnot("COND_CASE not found in condition column" =
              COND_CASE %in% HIV1@meta.data[[COND_COL]],
            "COND_CTRL not found in condition column" =
              COND_CTRL %in% HIV1@meta.data[[COND_COL]])
  
  # ----------------------------- 1. QC TABLES ----------------------------------
  HIV1@meta.data %>%
    dplyr::group_by(atac_clusters) %>%
    dplyr::summarise(n = dplyr::n(),
                     frac_dbl = mean(dbl_class == "doublet"),
                     rna_umi = median(nCount_RNA), genes = median(nFeature_RNA),
                     atac = median(nCount_ATAC), mt = median(percent.mt)) %>%
    tibble::as_tibble() %>% print(n = 30)
  
  print(table(HIV1$atac_clusters, HIV1@meta.data[[COND_COL]]))
  if (!is.null(DONOR_COL)) print(table(HIV1$atac_clusters, HIV@meta.data[[DONOR_COL]]))
  
  # ----------------------------- 2. FILTER (subset only, no re-clustering) -----
  DefaultAssay(HIV1) <- "RNA"
  Idents(HIV1) <- "atac_clusters"

  markers_9 <- FindMarkers(HIV1, ident.1 = "9", only.pos = TRUE, min.pct = 0.25)
  head(markers_9, 40)
  
  markers_18 <- FindMarkers(HIV1, ident.1 = "18", only.pos = TRUE, min.pct = 0.25)
  head(markers_18, 10)
  markers_10 <- FindMarkers(HIV1, ident.1 = "10", only.pos = TRUE, min.pct = 0.25)
  head(markers_10, 10)
  markers_20 <- FindMarkers(HIV1, ident.1 = "20", only.pos = TRUE, min.pct = 0.25)
  head(markers_20, 10)
  markers_14 <- FindMarkers(HIV1, ident.1 = "14", only.pos = TRUE, min.pct = 0.25)
  head(markers_14, 10)
  # ----------------------------- 3. ANNOTATION ---------------------------------
  # PROVISIONAL: edit after looking at the dot plots. "Unassigned" is excluded from DE.
  cluster_map <- c(
    "0"  = "CD4 T",       "1"  = "CD4 T",       "2"  = "CD8 T",       "3"  = "CD4 T",
    "4"  = "Cycling T",   "5"  = "CD4 T",       "6"  = "B cell",      "7"  = "CD8 T",
    "8"  = "DC",          "9"  = "Low_quality",  "10" = "cDC2",  "11" = "B cell",
    "12" = "DC",          "13" = "T innate-like","14" = "Monocyte", "15" = "CD4 T",
    "16" = "T innate-like","17" = "DC",         "18" = "Activated T", "19" = "Unassigned",
    "20" = "Activated T"
  )
  # Set identity to your new cell type column
  Idents(HIV_f) <- "celltype"
  
  # Generate the UMAP DimPlot grouped by cell type with labels
  DimPlot(HIV_f, reduction = "umap.atac", group.by = "celltype", label = TRUE, repel = TRUE) +
    ggtitle("ATAC UMAP - Cell Type Annotations") +
    theme_minimal()
  
  # Save the plot
  ggsave("results/umap_atac_annotated_celltypes.pdf", p_dim, width = 8, height = 7)
  print(p_dim)
  HIV_f<-HIV1
  missing_cl <- setdiff(levels(HIV_f$atac_clusters), names(cluster_map))
  if (length(missing_cl)) stop("Add to cluster_map: ", paste(missing_cl, collapse = ", "))
  
  HIV_f$celltype <- unname(cluster_map[as.character(HIV_f$atac_clusters)])
  print(table(HIV_f$atac_clusters, HIV_f$celltype))
  
  DefaultAssay(HIV_f) <- "RNA"
  p <- DotPlot(HIV_f, group.by = "atac_clusters",
               features = c("CD3E","CD4","CD8A","IL7R","CCR7","SELL","GZMK","GZMH",
                            "NKG7","GNLY","NCAM1","TRDC","MS4A1","CD79A","MKI67")) +
    RotatedAxis()
  ggsave("results/dotplot_lineage_markers.pdf", p, width = 9, height = 6)
  
  p <- DotPlot(HIV_f, group.by = "atac_clusters",
               features = c("CD1C","CLEC10A","FCER1A","CLEC9A","XCR1","IRF8",
                            "LAMP3","CCR7","CCL19","IDO1","IL3RA","LILRA4","TCF4",
                            "CD14","S100A8","HLA-DRA")) + RotatedAxis()
  ggsave("results/dotplot_DC_markers.pdf", p, width = 9, height = 6)
  
  DimPlot(HIV_f, reduction = "umap.atac", group.by = "celltype", label = TRUE)
  ggsave("results/umap_atac_celltype.pdf", width = 7, height = 6)
  
  # Every cell type needs cells from both conditions
  print(HIV_f@meta.data %>%
          dplyr::filter(celltype != "Unassigned") %>%
          dplyr::group_by(celltype, .data[[COND_COL]]) %>%
          dplyr::summarise(n_cells = dplyr::n(),
                           n_samples = dplyr::n_distinct(.data[[SAMPLE_COL]]),
                           .groups = "drop") %>%
          tibble::as_tibble(), n = 50)
  
  # ----------------------------- 4. PSEUDOBULK RNA DE --------------------------
  # Seurat v5: layer = "counts".  Seurat v4: use slot = "counts".
  counts <- GetAssayData(HIV_f, assay = "RNA", layer = "counts")
  
  md <- HIV_f@meta.data
  stopifnot(identical(colnames(counts), rownames(md)))
  md$key <- paste(md$celltype, md[[SAMPLE_COL]], sep = "||")
  keys <- factor(md$key)
  
  # cells x pseudobulks indicator matrix (no model.matrix needed)
  mm <- sparseMatrix(i = seq_along(keys), j = as.integer(keys), x = 1,
                     dims = c(length(keys), nlevels(keys)),
                     dimnames = list(colnames(counts), levels(keys)))
  pb <- counts %*% mm                         # genes x pseudobulks
  
  meta_pb <- md %>%
    dplyr::group_by(key) %>%
    dplyr::summarise(celltype  = dplyr::first(celltype),
                     sample    = dplyr::first(.data[[SAMPLE_COL]]),
                     condition = dplyr::first(.data[[COND_COL]]),
                     donor     = if (is.null(DONOR_COL)) NA_character_
                     else dplyr::first(.data[[DONOR_COL]]),
                     ncells    = dplyr::n(), .groups = "drop") %>%
    as.data.frame()
  rownames(meta_pb) <- meta_pb$key
  meta_pb <- meta_pb[colnames(pb), ]
  meta_pb$condition <- relevel(factor(meta_pb$condition), ref = COND_CTRL)
  
  # Replication check: pseudobulks per cell type x condition (>= MIN_CELLS)
  cat("\nPseudobulks per celltype x condition (>= MIN_CELLS cells):\n")
  print(table(meta_pb$celltype[meta_pb$ncells >= MIN_CELLS],
              meta_pb$condition[meta_pb$ncells >= MIN_CELLS]))
  cat("\nMedian cells per pseudobulk:\n")
  print(meta_pb %>% dplyr::group_by(celltype, condition) %>%
          dplyr::summarise(median_cells = median(ncells), .groups = "drop") %>%
          tibble::as_tibble(), n = 50)
  
  have_ashr <- requireNamespace("ashr", quietly = TRUE)
  if (!have_ashr) message("ashr not installed: lfc_shrunk will be NA (p-values unaffected).")
  
  run_pb_de <- function(ct) {
    cd <- meta_pb[meta_pb$celltype == ct & meta_pb$ncells >= MIN_CELLS, ]
    n_per <- table(cd$condition)
    if (!all(c(COND_CASE, COND_CTRL) %in% names(n_per)) ||
        any(n_per[c(COND_CASE, COND_CTRL)] < 2)) {
      message("Skipping ", ct, ": needs >= 2 pseudobulks per condition"); return(NULL)
    }
    m <- round(as.matrix(pb[, rownames(cd), drop = FALSE]))
    des <- if (!is.null(DONOR_COL) && !anyNA(cd$donor)) ~ donor + condition else ~ condition
    dds <- DESeqDataSetFromMatrix(m, cd, design = des)
    dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
    dds <- DESeq(dds, quiet = TRUE)
    
    res_raw <- results(dds, contrast = c("condition", COND_CASE, COND_CTRL))
    res <- as.data.frame(res_raw)                 # raw LFC: use for thresholds
    res$lfc_shrunk <- NA_real_                    # shrunk LFC: use for ranking/plots
    if (have_ashr) {
      res$lfc_shrunk <- lfcShrink(dds, contrast = c("condition", COND_CASE, COND_CTRL),
                                  res = res_raw, type = "ashr", quiet = TRUE)$log2FoldChange
    }
    res$gene <- rownames(res); res$celltype <- ct
    res
  }
  
  cts <- setdiff(unique(HIV_f$celltype), "Unassigned")
  de_list <- setNames(lapply(cts, run_pb_de), cts)
  de_all  <- dplyr::bind_rows(de_list)
  write.csv(de_all, "E:/results/pseudobulk_DE_by_celltype.csv", row.names = FALSE)
  
  cat("\nCell types skipped (not testable):",
      paste(names(de_list)[vapply(de_list, is.null, logical(1))], collapse = ", "), "\n")
  
  print(de_all %>% dplyr::group_by(celltype) %>%
          dplyr::summarise(n_up = sum(padj < 0.05 & log2FoldChange > 0, na.rm = TRUE),
                           n_down = sum(padj < 0.05 & log2FoldChange < 0, na.rm = TRUE)) %>%
          tibble::as_tibble())
  
  # Known ISGs, for a quick sanity check
  isg <- c("ISG15","IFI6","IFI44","IFI44L","IFI27","MX1","MX2","USP18","OAS1","IFIT1","IFIT3")
  print(de_all %>% dplyr::filter(gene %in% isg) %>%
          dplyr::select(celltype, gene, log2FoldChange, lfc_shrunk, padj) %>%
          dplyr::arrange(celltype, dplyr::desc(log2FoldChange)) %>%
          tibble::as_tibble(), n = 60)
  
  # ISG module score per pseudobulk: more power than single genes with few replicates
  isg_set <- intersect(c("ISG15","IFI6","IFI44","IFI44L","IFI27","MX1","MX2","USP18",
                         "OAS1","IFIT1","IFIT3","IFIT2","OAS2","OAS3","RSAD2","IFITM3","STAT1"),
                       rownames(pb))
  lib  <- Matrix::colSums(pb)
  cpm  <- t(t(as.matrix(pb)) / lib) * 1e6
  meta_pb$isg_score <- colMeans(log2(cpm[isg_set, rownames(meta_pb), drop = FALSE] + 1))
  ggplot(meta_pb[meta_pb$ncells >= MIN_CELLS & meta_pb$celltype != "Unassigned", ],
         aes(condition, isg_score, color = condition)) +
    geom_point(size = 3) + facet_wrap(~celltype, scales = "free_y") + theme_bw()
  ggsave("results/isg_score_by_celltype.pdf", width = 8, height = 6)
  
  # Optional cell-level exploratory DE (p-values inflated; ranking only)
  # Idents(HIV_f) <- "celltype"
  # sub <- subset(HIV_f, idents = "CD4 T"); Idents(sub) <- COND_COL
  # FindMarkers(sub, ident.1 = COND_CASE, ident.2 = COND_CTRL, test.use = "wilcox")
  
  # ----------------------------- 5. PRIMED-CHROMATIN TEST ----------------------
  DefaultAssay(HIV_f) <- "ATAC"
  annot <- Annotation(HIV_f)
  if (PROTEIN_CODING_ONLY && "gene_biotype" %in% colnames(mcols(annot))) {
    annot <- annot[annot$gene_biotype == "protein_coding"]
  }
  prom  <- promoters(annot, upstream = PROM_UP, downstream = PROM_DOWN)
  peaks <- granges(HIV_f[["ATAC"]])
  ov    <- findOverlaps(peaks, prom)
  peak2gene <- data.frame(peak = rownames(HIV_f[["ATAC"]])[queryHits(ov)],
                          gene = prom$gene_name[subjectHits(ov)]) %>%
    dplyr::distinct()
  
  # Funnel: where do genes get lost?
  funnel <- lapply(cts, function(ct) {
    de <- de_list[[ct]]
    if (is.null(de)) return(data.frame(celltype = ct, status = "skipped"))
    nonrp <- !grepl("^(RP[LS]|MRP[LS]|MT-)", de$gene)
    sel   <- !is.na(de$padj) & de$padj < PADJ_RNA & de$log2FoldChange > LFC_RNA &
      de$baseMean > MIN_BASEMEAN & nonrp
    data.frame(celltype = ct, status = "ok", n_tested = nrow(de),
               sig_up = sum(!is.na(de$padj) & de$padj < PADJ_RNA & de$log2FoldChange > 0),
               pass_filters = sum(sel),
               with_promoter_peak = sum(sel & de$gene %in% peak2gene$gene))
  })
  cat("\nFunnel:\n"); print(dplyr::bind_rows(funnel))
  
  primed_list <- list()
  for (ct in cts) {
    de <- de_list[[ct]]
    if (is.null(de)) next
    
    up <- de %>% dplyr::filter(!is.na(padj), padj < PADJ_RNA,
                               log2FoldChange > LFC_RNA, baseMean > MIN_BASEMEAN,
                               !grepl("^(RP[LS]|MRP[LS]|MT-)", gene))
    if (nrow(up) == 0) next
    
    ctrl <- rownames(meta_pb)[meta_pb$celltype == ct &
                                meta_pb$condition == COND_CTRL &
                                meta_pb$ncells >= MIN_CELLS]
    m <- as.matrix(pb[up$gene, ctrl, drop = FALSE])
    up$ctrl_cpm <- rowMeans(t(t(m) / lib[ctrl]) * 1e6)
    
    pk <- peak2gene %>% dplyr::filter(gene %in% up$gene)
    if (nrow(pk) == 0) next
    
    sub <- subset(HIV_f, subset = celltype == ct)
    Idents(sub) <- COND_COL
    # Cell-level LR test on promoter peaks only (fast; padj corrected over these peaks)
    da <- FindMarkers(sub, ident.1 = COND_CASE, ident.2 = COND_CTRL,
                      features = unique(pk$peak), test.use = "LR",
                      latent.vars = "nCount_ATAC",
                      min.pct = 0, logfc.threshold = 0, verbose = FALSE)
    da$peak <- rownames(da)
    
    prim <- pk %>%
      dplyr::left_join(da, by = "peak") %>%
      dplyr::mutate(open_both = pct.1 >= OPEN_PCT & pct.2 >= OPEN_PCT,
                    stable    = p_val_adj > 0.05 & abs(avg_log2FC) < STABLE_LFC,
                    opens     = (p_val_adj < 0.05 & avg_log2FC > 0) |
                      (pct.1 >= OPEN_PCT & pct.2 < OPEN_PCT)) %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(primed = any(open_both & stable, na.rm = TRUE),
                       opens_with_stimulus = any(opens, na.rm = TRUE),
                       .groups = "drop") %>%
      dplyr::left_join(up %>% dplyr::select(gene, log2FoldChange, lfc_shrunk, padj, ctrl_cpm),
                       by = "gene") %>%
      dplyr::mutate(category = dplyr::case_when(
        primed & ctrl_cpm < CPM_SILENT ~ "primed (silent at baseline)",
        primed                         ~ "open, already expressed",
        opens_with_stimulus            ~ "chromatin opens with stimulus",
        TRUE                           ~ "other / weak peak"),
        celltype = ct)
    primed_list[[ct]] <- prim
  }
  
  primed <- dplyr::bind_rows(primed_list)
  if (nrow(primed) > 0) {
    print(table(primed$celltype, primed$category))
    write.csv(primed, "results/primed_genes.csv", row.names = FALSE)
    print(primed %>% dplyr::arrange(celltype, dplyr::desc(log2FoldChange)) %>%
            dplyr::select(celltype, gene, log2FoldChange, ctrl_cpm, category) %>%
            tibble::as_tibble(), n = 60)
  } else {
    message("No genes passed the filters. Relax LFC_RNA / PADJ_RNA or check replication above.")
  }
  
  # Coverage plots for top induced genes in the best-powered cell type
  if (nrow(primed) > 0) {
    ct_show <- names(sort(table(primed$celltype), decreasing = TRUE))[1]
    top_genes <- primed %>% dplyr::filter(celltype == ct_show) %>%
      dplyr::arrange(padj) %>% head(4) %>% dplyr::pull(gene)
    sub <- subset(HIV_f, subset = celltype == ct_show)
    Idents(sub) <- COND_COL
    for (g in top_genes) {
      p <- CoveragePlot(sub, region = g, features = g, expression.assay = "RNA",
                        extend.upstream = 5000, extend.downstream = 2000)
      ggsave(sprintf("results/coverage_%s_%s.pdf", gsub("[ /]", "_", ct_show), g),
             p, width = 8, height = 6)
    }
  }
  annot <- Annotation(HIV_f)
  table(annot$type)
  if (PROTEIN_CODING_ONLY && "gene_biotype" %in% colnames(mcols(annot))) {
    annot <- annot[annot$gene_biotype == "protein_coding"]
  }
  gene_gr <- unlist(range(split(annot, annot$gene_name)))   # one range per gene, strand kept
  gene_gr$gene_name <- names(gene_gr)
  prom <- promoters(gene_gr, upstream = PROM_UP, downstream = PROM_DOWN)
  prom
  # ----------------------------- 6. SAVE ---------------------------------------
  saveRDS(HIV_f, "results/HIV_filtered_annotated.rds")
  saveRDS(list(pb = pb, meta_pb = meta_pb, de_list = de_list),
          "results/pseudobulk_objects.rds")
  
  
  
  dir.create("E:/bed", showWarnings = FALSE)
  DefaultAssay(HIV_f) <- "ATAC"
  counts <- GetAssayData(HIV_f, assay = "ATAC", layer = "counts")   # Seurat v4: slot = "counts"
  pk <- granges(HIV_f[["ATAC"]])
  OPEN_PCT <- 0.05          # keep peaks detected in >= 5% of cells of that condition
  
  for (cond in c("HIV", "Mock")) {
    cells <- colnames(HIV_f)[HIV_f@meta.data[[COND_COL]] == cond]
    m     <- counts[, cells, drop = FALSE]
    pct   <- Matrix::rowMeans(m > 0)          # fraction of cells with the peak detected
    total <- Matrix::rowSums(m)               # total fragment counts in the peak
    keep  <- pct >= OPEN_PCT
    
    bed <- data.frame(chr   = as.character(seqnames(pk)),
                      start = start(pk) - 1,  # BED is 0-based
                      end   = end(pk),
                      name  = rownames(HIV_f[["ATAC"]]),
                      score = round(total),
                      strand = ".",
                      pct_cells = round(pct, 4))[keep, ]
    write.table(bed, sprintf("E:/bed/%s_peaks.bed", cond), sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    cat(cond, ":", nrow(bed), "peaks,", length(cells), "cells\n")
  }
  
  table(HIV_f@meta.data[[COND_COL]])
  HIV_f@meta.data %>% dplyr::group_by(.data[[COND_COL]]) %>%
    dplyr::summarise(cells = dplyr::n(), median_atac = median(nCount_ATAC)) %>%
    tibble::as_tibble()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
   
   
  library(Signac)
  library(GenomicRanges)
  
  # Convert peak names to genomic ranges
  peak_ranges <- StringToGRanges(
    rownames(DA_HIV_vs_Mock_cluster0),
    sep = c("-", "-")
  )
  
  # Keep significant peaks
  sig_peak_ranges <- peak_ranges[
    rownames(peak_ranges) %in% sig_peaks
  ]
  
  sig_peak_ranges
  
  
  
  
  
  
  
  
  
  
  # Set identity to cluster IDs
  Idents(cluster0_HIV_Mock) <- "seurat_clusters" # Adjust to your cluster metadata column name
  clusters <- levels(Idents(cluster0_HIV_Mock))  # Assuming 18 clusters (0 to 17)
  
  # Create an output directory for saving results
  dir.create("DA_Peaks_HIV_vs_Mock", showWarnings = FALSE)
  
  # Loop over all clusters
  for (cluster in clusters) {
    
    message(paste0("Processing Cluster ", cluster, "..."))
    
    # 1. Subset object to the specific cluster
    sub_obj <- subset(cluster0_HIV_Mock, idents = cluster)
    
    # Set identities within subset to the condition variable
    Idents(sub_obj) <- "condition" # Adjust to your HIV/Mock column name (e.g., 'group')
    
    # 2. Check if both HIV and Mock conditions exist in this cluster
    if (all(c("HIV", "Mock") %in% Idents(sub_obj))) {
      
      # 3. Find differentially accessible peaks (HIV vs Mock)
      da_peaks <- FindMarkers(
        object = sub_obj,
        ident.1 = "HIV",
        ident.2 = "Mock",
        only.pos = FALSE,             # Gets both up- and down-regulated peaks
        min.pct = 0.05,
        test.use = "LR",              # Logistic Regression (standard for ATAC)
        latent.vars = "nCount_ATAC"   # Adjusts for depth differences
      )
      
      # Add peak coordinates as a distinct column
      da_peaks$peak <- rownames(da_peaks)
      
      # 4. Save results to CSV and RDS
      write.csv(
        da_peaks, 
        file = paste0("DA_Peaks_HIV_vs_Mock/Cluster_", cluster, "_DA_peaks.csv"), 
        row.names = FALSE
      )
      
      saveRDS(
        da_peaks, 
        file = paste0("DA_Peaks_HIV_vs_Mock/Cluster_", cluster, "_DA_peaks.rds")
      )
      
      message(paste0("Saved results for Cluster ", cluster))
      
    } else {
      warning(paste0("Cluster ", cluster, " does not contain both HIV and Mock cells. Skipping."))
    }
  }
  
  
  table(HIV1$condition)
  Fragments(HIV1)
  
  Fragments(HIV1[["ATAC"]])
  frags <- Fragments(HIV1[["ATAC"]])[[1]]
  
  frags 
  
  
  head(frags@cells)
  head(colnames(HIV1))
  
  barcode <- sub(".*_", "", colnames(HIV1))
  
  length(frags@cells)
  sum(sub(".*_", "", colnames(HIV1)) %in% frags@cells)
  