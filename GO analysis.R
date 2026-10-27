# ============================================================
# HIV vs ssDNA (mock), RNA, all clusters
# Step 1: differential expression per cluster (HIV vs ssDNA)
# Step 2: GO (ORA) + GSEA per cluster, human
# Step 3: summary tables and overview plots across clusters
# PNG output (300 dpi), large fonts, wrapped term labels.
# ============================================================
library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
library(clusterProfiler); library(enrichplot); library(org.Hs.eg.db)

# ---- Settings -------------------------------------------------
obj_path      <- "E:/multiome_with_Marie_celltype_atac.rds"
rna_assay     <- "RNA"
cluster_col   <- "seurat_clusters"
condition_col <- "condition"
hiv_label     <- "ssDNA"
mock_label    <- "Mock"

out_root <- "E:/ssDNA_vs_Mock_GO"          # keep this path SHORT (Windows 260-char limit)
de_dir   <- file.path(out_root, "DE")
gdir     <- file.path(out_root, "GO")

min_cells_per_condition <- 10
orgdb <- org.Hs.eg.db
ont <- "BP"; padj_cut <- 0.05; lfc_cut <- 0.25
fs <- 20                                    # base font size for all plots
wrap_n <- 40                                # characters per line for GO term labels

lab_hiv  <- paste("Higher in", hiv_label)
lab_mock <- paste("Higher in", mock_label)

for (d in c(out_root, de_dir, gdir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(gdir)) stop("Cannot create output folder: ", gdir)

safe <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

big <- theme(
  axis.text    = element_text(size = fs),
  axis.title   = element_text(size = fs + 5),
  plot.title   = element_text(size = fs + 5, face = "bold"),
  legend.text  = element_text(size = fs - 2),
  legend.title = element_text(size = fs),
  strip.text   = element_text(size = fs)
)
wrapy <- scale_y_discrete(labels = scales::label_wrap(wrap_n))

# ============================================================
# STEP 1: differential expression per cluster
# ============================================================
run_de_all <- function() {
  obj <- readRDS(obj_path)
  DefaultAssay(obj) <- rna_assay
  obj <- tryCatch(JoinLayers(obj, assay = rna_assay), error = function(e) obj)
  if (!"data" %in% Layers(obj, assay = rna_assay)) obj <- NormalizeData(obj)
  
  md <- obj[[]]
  clusters <- sort(unique(as.character(md[[cluster_col]])))
  message("Found clusters: ", paste(clusters, collapse = ", "))
  
  for (cl in clusters) {
    out_csv <- file.path(de_dir, paste0("Cluster_", safe(cl), "_ssDNA_vs_Mock.csv"))
    if (file.exists(out_csv)) { message("skip DE (exists): cluster ", cl); next }
    
    cells <- rownames(md)[as.character(md[[cluster_col]]) == cl &
                            as.character(md[[condition_col]]) %in% c(hiv_label, mock_label)]
    cnt <- table(as.character(md[cells, condition_col]))
    if (length(cnt) < 2 || any(cnt < min_cells_per_condition)) {
      message("Skipping cluster ", cl, ": too few cells (", paste(names(cnt), cnt, collapse = ", "), ")")
      next
    }
    message("DE cluster ", cl, ": ", paste(names(cnt), cnt, collapse = ", "))
    
    sub <- subset(obj, cells = cells)
    Idents(sub) <- condition_col
    res <- tryCatch(
      FindMarkers(sub, assay = rna_assay, ident.1 = hiv_label, ident.2 = mock_label,
                  test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0),
      error = function(e) { message("FindMarkers failed: ", e$message); NULL })
    if (is.null(res) || nrow(res) == 0) next
    
    res$gene <- rownames(res)
    write.csv(res, out_csv, row.names = FALSE)   # positive avg_log2FC = higher in HIV
    rm(sub); gc()
  }
  rm(obj); gc()
}

# ============================================================
# STEP 2: GO (ORA) + GSEA for one DE file
# ============================================================
run_go_gsea <- function(csv_name) {
  tag <- safe(sub("\\.csv$", "", csv_name))
  out <- file.path(gdir, tag); dir.create(out, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(out)) stop("Cannot create folder (path too long?): ", out)
  
  sv <- function(pl, nm, w, h)
    ggsave(file.path(out, nm), pl, width = w, height = h, units = "in", dpi = 300, bg = "white")
  
  deg <- read.csv(file.path(de_dir, csv_name), check.names = FALSE)
  
  sig <- deg %>%
    filter(p_val_adj < padj_cut, abs(avg_log2FC) >= lfc_cut) %>%
    mutate(direction = ifelse(avg_log2FC > 0, lab_hiv, lab_mock)) %>%
    arrange(desc(avg_log2FC))
  write.csv(sig, file.path(out, "significant_DEGs.csv"), row.names = FALSE)
  cat("\n==", csv_name, "==\n"); print(table(sig$direction))
  
  # GO enrichment (ORA)
  gs <- list(sig$gene[sig$direction == lab_hiv], sig$gene[sig$direction == lab_mock])
  names(gs) <- c(lab_hiv, lab_mock)
  gs <- gs[lengths(gs) >= 5]
  ego <- if (length(gs)) tryCatch(
    compareCluster(gs, fun = "enrichGO", OrgDb = orgdb, keyType = "SYMBOL", ont = ont,
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05),
    error = function(e) { message("ORA error: ", e$message); NULL })
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego <- tryCatch(simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min),
                    error = function(e) ego)
    res <- as.data.frame(ego)
    write.csv(res, file.path(out, "GO_ORA_terms.csv"), row.names = FALSE)
    
    sv(dotplot(ego, showCategory = 10, font.size = fs) + big + wrapy +
         ggtitle(paste("GO", ont)), "GO_ORA_dotplot.png", 16, 15)
    
    for (d in unique(res$Cluster)) {
      x <- res %>% filter(Cluster == d) %>% slice_min(p.adjust, n = 15, with_ties = FALSE) %>%
        mutate(Description = factor(Description, levels = rev(Description)))
      sv(ggplot(x, aes(Count, Description, fill = p.adjust)) + geom_col() +
           theme_bw(base_size = fs) + big + wrapy +
           scale_fill_gradient(low = "red", high = "blue") +
           labs(title = paste("GO", ont, "-", d), y = NULL),
         paste0("GO_ORA_barplot_", safe(d), ".png"), 16, 12)
    }
  } else message("No significant GO terms (ORA)")
  
  # GSEA, ranked by avg_log2FC (positive = higher in ssDNA)
  rk <- setNames(deg$avg_log2FC, deg$gene)
  rk <- sort(rk[!is.na(rk) & !duplicated(names(rk))], decreasing = TRUE)
  gsea <- tryCatch(
    gseGO(rk, OrgDb = orgdb, keyType = "SYMBOL", ont = ont, minGSSize = 10, maxGSSize = 500,
          pvalueCutoff = 0.05, eps = 0, verbose = FALSE),
    error = function(e) { message("GSEA error: ", e$message); NULL })
  
  if (!is.null(gsea) && nrow(gsea@result) > 0) {
    write.csv(gsea@result, file.path(out, "GSEA_GO_results.csv"), row.names = FALSE)
    
    sv(dotplot(gsea, showCategory = 10, split = ".sign", font.size = fs) + big + wrapy +
         facet_grid(. ~ .sign) + ggtitle(paste("GSEA", ont, "(activated = higher in", hiv_label, ")")),
       "GSEA_dotplot.png", 22, 15)
    
    sv(ridgeplot(gsea, showCategory = 12) + big + wrapy, "GSEA_ridgeplot.png", 16, 14)
    
    r   <- gsea@result %>% arrange(desc(NES))
    ids <- match(c(head(r$ID[r$NES > 0], 3), tail(r$ID[r$NES < 0], 3)), gsea@result$ID)
    ids <- ids[!is.na(ids)]
    if (length(ids) > 0) {
      g <- gsea
      d <- g@result$Description
      d <- ifelse(nchar(d) > 45, paste0(substr(d, 1, 42), "..."), d)
      g@result$Description <- paste0(d, "  (NES=", round(g@result$NES, 2),
                                     ", padj=", signif(g@result$p.adjust, 2), ")")
      sv(gseaplot2(g, geneSetID = ids, pvalue_table = FALSE, base_size = fs + 2),
         "GSEA_curves.png", 22, 14)
    }
  } else message("No significant GSEA terms")
  
  message("Done: ", out)
  invisible(NULL)
}

# ============================================================
# STEP 3: summary across clusters
# ============================================================
summarize_all <- function() {
  pdir <- file.path(gdir, "summary_plots"); dir.create(pdir, showWarnings = FALSE)
  dirs <- setdiff(list.dirs(gdir, recursive = FALSE), pdir)
  sv2  <- function(pl, nm, w, h)
    ggsave(file.path(pdir, nm), pl, width = w, height = h, units = "in", dpi = 300,
           bg = "white", limitsize = FALSE)
  
  read_all <- function(fn) bind_rows(lapply(dirs, function(d) {
    f <- file.path(d, fn); if (!file.exists(f)) return(NULL)
    x <- read.csv(f, check.names = FALSE)
    x$cluster_id <- sub("_ssDNA_vs_Mock$", "", sub("^Cluster_", "", basename(d)))
    x }))
  
  DEG <- read_all("significant_DEGs.csv")
  GO  <- read_all("GO_ORA_terms.csv")
  GS  <- read_all("GSEA_GO_results.csv")
  if (!is.null(DEG) && nrow(DEG)) write.csv(DEG, file.path(out_root, "ALL_significant_DEGs.csv"), row.names = FALSE)
  if (!is.null(GO)  && nrow(GO))  write.csv(GO,  file.path(out_root, "ALL_GO_ORA_terms.csv"),     row.names = FALSE)
  if (!is.null(GS)  && nrow(GS))  write.csv(GS,  file.path(out_root, "ALL_GSEA_GO_results.csv"),  row.names = FALSE)
  
  dir_cols <- setNames(c("firebrick", "steelblue"), c(lab_hiv, lab_mock))
  
  # DEG counts per cluster
  if (!is.null(DEG) && nrow(DEG)) {
    sv2(ggplot(count(DEG, cluster_id, direction), aes(cluster_id, n, fill = direction)) +
          geom_col(position = "dodge") + theme_bw(base_size = fs) + big +
          scale_fill_manual(values = dir_cols) +
          labs(x = "Cluster", y = "Significant genes", fill = NULL) +
          theme(legend.position = "bottom"),
        "DEG_counts.png", 14, 8)
  }
  
  # GO dotplot across clusters, one per direction
  if (!is.null(GO) && nrow(GO)) {
    for (s in unique(GO$Cluster)) {
      d <- GO %>% filter(Cluster == s)
      top <- d %>% group_by(cluster_id) %>% slice_min(p.adjust, n = 3, with_ties = FALSE) %>%
        pull(Description) %>% unique()
      d <- d %>% filter(Description %in% top)
      sv2(ggplot(d, aes(cluster_id, Description, size = Count, color = p.adjust)) + geom_point() +
            scale_color_gradient(low = "red", high = "blue") +
            theme_bw(base_size = fs) + big + wrapy +
            labs(title = paste("GO", ont, "-", s), x = "Cluster", y = NULL),
          paste0("GO_dotplot_", safe(s), ".png"), 16, max(9, 0.5 * length(top) + 5))
    }
  }
  message("Summary saved in: ", pdir)
}

# ============================================================
# Run everything
# ============================================================
run_de_all()

files <- list.files(de_dir, pattern = "^Cluster_.*\\.csv$")
cat("DE files to process:", length(files), "\n")

failed <- character(0)
t0 <- Sys.time()
for (i in seq_along(files)) {
  f    <- files[i]
  tag  <- safe(sub("\\.csv$", "", f))
  mark <- file.path(gdir, tag, "done.txt")
  if (file.exists(mark)) { message(sprintf("[%d/%d] skip (already done): %s", i, length(files), f)); next }
  
  message(sprintf("\n[%d/%d] %s   (elapsed: %.1f min)", i, length(files), f,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  ok <- tryCatch({ run_go_gsea(f); TRUE }, error = function(e) {
    message("FAILED: ", f, " -> ", e$message); FALSE })
  if (ok) file.create(mark) else failed <- c(failed, f)
  graphics.off(); gc()
}

cat("\nFinished. Total time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
cat("Failed files:", length(failed), "\n")
if (length(failed)) {
  print(failed)
  write.csv(data.frame(file = failed), file.path(out_root, "failed_files.csv"), row.names = FALSE)
}

summarize_all()
message("All done. Results in: ", out_root)
































# ============================================================
# HIV vs Mock, RNA, all clusters
# Step 1: differential expression per cluster (HIV vs Mock)
# Step 2: GO (ORA) + GSEA per cluster, human
#         (ribosomal and mitochondrial genes removed before GO/GSEA)
# Step 3: summary tables and overview plots across clusters
# PNG output (300 dpi), large fonts, wrapped term labels.
# ============================================================
library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
library(clusterProfiler); library(enrichplot); library(org.Hs.eg.db)

# ---- Settings -------------------------------------------------
obj_path      <- "E:/multiome_with_Marie_celltype_atac.rds"
rna_assay     <- "RNA"
cluster_col   <- "seurat_clusters"
condition_col <- "condition"
hiv_label     <- "HIV"
mock_label    <- "ssDNA"

out_root <- "E:/HIV_vs_ssDNA_GO"          # keep this path SHORT (Windows 260-char limit)
de_dir   <- file.path(out_root, "DEG/DEG")
gdir     <- file.path(out_root, "GO")

# TRUE : run the DE step (loads the Seurat object; clusters with an existing DE file are skipped)
# FALSE: use the stored DE files in de_dir as they are (the object is NOT loaded, much faster)
run_de <- FALSE

min_cells_per_condition <- 10
orgdb <- org.Hs.eg.db
ont <- "BP"; padj_cut <- 0.05; lfc_cut <- 0.25
fs <- 20                                    # base font size for all plots
wrap_n <- 40                                # characters per line for GO term labels

# Genes removed before GO and GSEA (the saved DE files are NOT changed):
#   MT-*            mitochondrial-encoded genes
#   RPL*/RPS*       cytosolic ribosomal proteins (incl. RPLP0/1/2, RPSA, ribosomal pseudogenes)
#   MRPL*/MRPS*     mitochondrial ribosomal proteins (remove this line if you want to keep them)
drop_pattern <- "^MT-|^RP[SL][0-9]|^RPLP[0-2]$|^RPSA$|^MRP[SL][0-9]"

lab_hiv  <- paste("Higher in", hiv_label)
lab_mock <- paste("Higher in", mock_label)

for (d in c(out_root, de_dir, gdir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(gdir)) stop("Cannot create output folder: ", gdir)

safe <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

big <- theme(
  axis.text    = element_text(size = fs),
  axis.title   = element_text(size = fs + 5),
  plot.title   = element_text(size = fs + 5, face = "bold"),
  legend.text  = element_text(size = fs - 2),
  legend.title = element_text(size = fs),
  strip.text   = element_text(size = fs)
)
wrapy <- scale_y_discrete(labels = scales::label_wrap(wrap_n))

# ============================================================
# STEP 1: differential expression per cluster
# ============================================================
run_de_all <- function() {
  obj <- readRDS(obj_path)
  DefaultAssay(obj) <- rna_assay
  obj <- tryCatch(JoinLayers(obj, assay = rna_assay), error = function(e) obj)
  if (!"data" %in% Layers(obj, assay = rna_assay)) obj <- NormalizeData(obj)
  
  md <- obj[[]]
  clusters <- sort(unique(as.character(md[[cluster_col]])))
  message("Found clusters: ", paste(clusters, collapse = ", "))
  
  for (cl in clusters) {
    out_csv <- file.path(de_dir, paste0("Cluster_", safe(cl), "_HIV_vs_Mock.csv"))
    if (file.exists(out_csv)) { message("skip DE (exists): cluster ", cl); next }
    
    cells <- rownames(md)[as.character(md[[cluster_col]]) == cl &
                            as.character(md[[condition_col]]) %in% c(hiv_label, mock_label)]
    cnt <- table(as.character(md[cells, condition_col]))
    if (length(cnt) < 2 || any(cnt < min_cells_per_condition)) {
      message("Skipping cluster ", cl, ": too few cells (", paste(names(cnt), cnt, collapse = ", "), ")")
      next
    }
    message("DE cluster ", cl, ": ", paste(names(cnt), cnt, collapse = ", "))
    
    sub <- subset(obj, cells = cells)
    Idents(sub) <- condition_col
    res <- tryCatch(
      FindMarkers(sub, assay = rna_assay, ident.1 = hiv_label, ident.2 = mock_label,
                  test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0),
      error = function(e) { message("FindMarkers failed: ", e$message); NULL })
    if (is.null(res) || nrow(res) == 0) next
    
    res$gene <- rownames(res)
    write.csv(res, out_csv, row.names = FALSE)   # positive avg_log2FC = higher in HIV
    rm(sub); gc()
  }
  rm(obj); gc()
}

# ============================================================
# STEP 2: GO (ORA) + GSEA for one DE file
# ============================================================
run_go_gsea <- function(csv_name) {
  tag <- safe(sub("\\.csv$", "", csv_name))
  out <- file.path(gdir, tag); dir.create(out, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(out)) stop("Cannot create folder (path too long?): ", out)
  
  sv <- function(pl, nm, w, h)
    ggsave(file.path(out, nm), pl, width = w, height = h, units = "in", dpi = 300, bg = "white")
  
  deg <- read.csv(file.path(de_dir, csv_name), check.names = FALSE)
  
  # remove ribosomal / mitochondrial genes before ORA and GSEA
  n0  <- nrow(deg)
  deg <- deg[!grepl(drop_pattern, deg$gene, ignore.case = TRUE), ]
  message("Removed ", n0 - nrow(deg), " ribosomal/mitochondrial genes; ", nrow(deg), " remain")
  
  sig <- deg %>%
    filter(p_val_adj < padj_cut, abs(avg_log2FC) >= lfc_cut) %>%
    mutate(direction = ifelse(avg_log2FC > 0, lab_hiv, lab_mock)) %>%
    arrange(desc(avg_log2FC))
  write.csv(sig, file.path(out, "significant_DEGs.csv"), row.names = FALSE)
  cat("\n==", csv_name, "==\n"); print(table(sig$direction))
  
  # GO enrichment (ORA)
  gs <- list(sig$gene[sig$direction == lab_hiv], sig$gene[sig$direction == lab_mock])
  names(gs) <- c(lab_hiv, lab_mock)
  gs <- gs[lengths(gs) >= 5]
  ego <- if (length(gs)) tryCatch(
    compareCluster(gs, fun = "enrichGO", OrgDb = orgdb, keyType = "SYMBOL", ont = ont,
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05),
    error = function(e) { message("ORA error: ", e$message); NULL })
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego <- tryCatch(simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min),
                    error = function(e) ego)
    res <- as.data.frame(ego)
    write.csv(res, file.path(out, "GO_ORA_terms.csv"), row.names = FALSE)
    
    sv(dotplot(ego, showCategory = 10, font.size = fs) + big + wrapy +
         ggtitle(paste("GO", ont)), "GO_ORA_dotplot.png", 16, 15)
    
    for (d in unique(res$Cluster)) {
      x <- res %>% filter(Cluster == d) %>% slice_min(p.adjust, n = 15, with_ties = FALSE) %>%
        mutate(Description = factor(Description, levels = rev(Description)))
      sv(ggplot(x, aes(Count, Description, fill = p.adjust)) + geom_col() +
           theme_bw(base_size = fs) + big + wrapy +
           scale_fill_gradient(low = "red", high = "blue") +
           labs(title = paste("GO", ont, "-", d), y = NULL),
         paste0("GO_ORA_barplot_", safe(d), ".png"), 16, 12)
    }
  } else message("No significant GO terms (ORA)")
  
  # GSEA, ranked by avg_log2FC (positive = higher in HIV)
  rk <- setNames(deg$avg_log2FC, deg$gene)
  rk <- sort(rk[!is.na(rk) & !duplicated(names(rk))], decreasing = TRUE)
  gsea <- tryCatch(
    gseGO(rk, OrgDb = orgdb, keyType = "SYMBOL", ont = ont, minGSSize = 10, maxGSSize = 500,
          pvalueCutoff = 0.05, eps = 0, verbose = FALSE),
    error = function(e) { message("GSEA error: ", e$message); NULL })
  
  if (!is.null(gsea) && nrow(gsea@result) > 0) {
    write.csv(gsea@result, file.path(out, "GSEA_GO_results.csv"), row.names = FALSE)
    
    sv(dotplot(gsea, showCategory = 10, split = ".sign", font.size = fs) + big + wrapy +
         facet_grid(. ~ .sign) + ggtitle(paste("GSEA", ont, "(activated = higher in", hiv_label, ")")),
       "GSEA_dotplot.png", 22, 15)
    
    sv(ridgeplot(gsea, showCategory = 12) + big + wrapy, "GSEA_ridgeplot.png", 16, 14)
    
    r   <- gsea@result %>% arrange(desc(NES))
    ids <- match(c(head(r$ID[r$NES > 0], 3), tail(r$ID[r$NES < 0], 3)), gsea@result$ID)
    ids <- ids[!is.na(ids)]
    if (length(ids) > 0) {
      g <- gsea
      d <- g@result$Description
      d <- ifelse(nchar(d) > 45, paste0(substr(d, 1, 42), "..."), d)
      g@result$Description <- paste0(d, "  (NES=", round(g@result$NES, 2),
                                     ", padj=", signif(g@result$p.adjust, 2), ")")
      sv(gseaplot2(g, geneSetID = ids, pvalue_table = FALSE, base_size = fs + 2),
         "GSEA_curves.png", 22, 14)
    }
  } else message("No significant GSEA terms")
  
  message("Done: ", out)
  invisible(NULL)
}

# ============================================================
# STEP 3: summary across clusters
# ============================================================
summarize_all <- function() {
  pdir <- file.path(gdir, "summary_plots"); dir.create(pdir, showWarnings = FALSE)
  dirs <- setdiff(list.dirs(gdir, recursive = FALSE), pdir)
  sv2  <- function(pl, nm, w, h)
    ggsave(file.path(pdir, nm), pl, width = w, height = h, units = "in", dpi = 300,
           bg = "white", limitsize = FALSE)
  
  read_all <- function(fn) bind_rows(lapply(dirs, function(d) {
    f <- file.path(d, fn); if (!file.exists(f)) return(NULL)
    x <- read.csv(f, check.names = FALSE)
    x$cluster_id <- sub(paste0("_HIV_vs_", mock_label, "$"), "", sub("^Cluster_", "", basename(d)))
    x }))
  
  DEG <- read_all("significant_DEGs.csv")
  GO  <- read_all("GO_ORA_terms.csv")
  GS  <- read_all("GSEA_GO_results.csv")
  if (!is.null(DEG) && nrow(DEG)) write.csv(DEG, file.path(out_root, "ALL_significant_DEGs.csv"), row.names = FALSE)
  if (!is.null(GO)  && nrow(GO))  write.csv(GO,  file.path(out_root, "ALL_GO_ORA_terms.csv"),     row.names = FALSE)
  if (!is.null(GS)  && nrow(GS))  write.csv(GS,  file.path(out_root, "ALL_GSEA_GO_results.csv"),  row.names = FALSE)
  
  dir_cols <- setNames(c("firebrick", "steelblue"), c(lab_hiv, lab_mock))
  
  # DEG counts per cluster
  if (!is.null(DEG) && nrow(DEG)) {
    sv2(ggplot(count(DEG, cluster_id, direction), aes(cluster_id, n, fill = direction)) +
          geom_col(position = "dodge") + theme_bw(base_size = fs) + big +
          scale_fill_manual(values = dir_cols) +
          labs(x = "Cluster", y = "Significant genes", fill = NULL) +
          theme(legend.position = "bottom"),
        "DEG_counts.png", 14, 8)
  }
  
  # GO dotplot across clusters, one per direction
  if (!is.null(GO) && nrow(GO)) {
    for (s in unique(GO$Cluster)) {
      d <- GO %>% filter(Cluster == s)
      top <- d %>% group_by(cluster_id) %>% slice_min(p.adjust, n = 3, with_ties = FALSE) %>%
        pull(Description) %>% unique()
      d <- d %>% filter(Description %in% top)
      sv2(ggplot(d, aes(cluster_id, Description, size = Count, color = p.adjust)) + geom_point() +
            scale_color_gradient(low = "red", high = "blue") +
            theme_bw(base_size = fs) + big + wrapy +
            labs(title = paste("GO", ont, "-", s), x = "Cluster", y = NULL),
          paste0("GO_dotplot_", safe(s), ".png"), 16, max(9, 0.5 * length(top) + 5))
    }
  }
  message("Summary saved in: ", pdir)
}

# ============================================================
# Run everything
# ============================================================
if (run_de) run_de_all() else message("Using stored DE files in: ", de_dir)

files <- list.files(de_dir, pattern = "^Cluster_.*\\.csv$")
cat("DE files to process:", length(files), "\n")
if (length(files) == 0) stop("No DE files found in ", de_dir,
                             ". Check de_dir, or set run_de <- TRUE.")

failed <- character(0)
t0 <- Sys.time()
for (i in seq_along(files)) {
  f    <- files[i]
  tag  <- safe(sub("\\.csv$", "", f))
  mark <- file.path(gdir, tag, "done.txt")
  if (file.exists(mark)) { message(sprintf("[%d/%d] skip (already done): %s", i, length(files), f)); next }
  
  message(sprintf("\n[%d/%d] %s   (elapsed: %.1f min)", i, length(files), f,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  ok <- tryCatch({ run_go_gsea(f); TRUE }, error = function(e) {
    message("FAILED: ", f, " -> ", e$message); FALSE })
  if (ok) file.create(mark) else failed <- c(failed, f)
  graphics.off(); gc()
}

cat("\nFinished. Total time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
cat("Failed files:", length(failed), "\n")
if (length(failed)) {
  print(failed)
  write.csv(data.frame(file = failed), file.path(out_root, "failed_files.csv"), row.names = FALSE)
}

summarize_all()
message("All done. Results in: ", out_root)